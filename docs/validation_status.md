# Validation Status

Last updated: November 12, 2026.

ARIA uses explicit validation boundaries so users can distinguish mature
workflows from beta analysis paths and implementation scaffolds.

## Authoritative modality tiers

This table MIRRORS the single source of truth — `MODALITY_VALIDATION` in
`aria/agents/orchestrator_agent.py` — and is verified against it by
`tests/test_docs_drift_guard.py` (the Docs Drift Guard). Edit the orchestrator
first; this table must match or the guard fails.

<!-- MODALITY_TIERS_TABLE_START -->
| Modality | Tier |
|---|---|
| scRNA | production |
| bulk_RNA | production |
| bulk_RNA_raw | production |
| scATAC | production |
| bulk_ATAC | production |
| ChIP | scaffold |
| CUT_AND_RUN | scaffold |
| CUT_AND_TAG | scaffold |
| HiC | scaffold |
<!-- MODALITY_TIERS_TABLE_END -->

**What "validated" means here.** ARIA's validated paths have been exercised on
controlled synthetic data and small real datasets, with reviewed reports and a
numerical ground-truth benchmark. This is **not** a claim of publication-grade
results for any specific study: a domain expert must still review the design, the
fitted model, and the conclusions before publication. ARIA is built to make that
review fast and honest — every claim is evidence-tiered, every analysis is
provenance-stamped, and missing or low-power results stay visible.

## Validated (controlled + small real datasets)

| Area | Status | Evidence |
|---|---|---|
| Bulk RNA-seq count matrix | Validated | Synthetic ground-truth DE benchmark (recall 1.0 / empirical FDR 0.0) + H9 three-condition workflow |
| Bulk RNA-seq raw FASTQ (`bulk_RNA_raw`, ADR-075) | Validated (`production`) | Universal multi-platform raw FASTQ engine validated end-to-end across `full_length` (SE/PE/multi-lane STAR 2-pass + `pysam` fail-closed auto-strand + `featureCounts`), `tag_3prime` (QuantSeq/SMART-3Seq poly(X) trimming + CPM export / TPM suppression), `tag_5prime` (CAGE/RAMPAGE), `umi_bulk` (`fastp --umi` + directional 1-Hamming `pysam` UMI collapse), and `pooled_plate_bulk` (BRB-seq/DRUG-seq/Prime-seq `STARsolo CB_UMI_Simple 1MM_Directional_UMItools`), plus H9 real chr21 raw FASTQ E2E (`c1_h9_fastq_e2e`), hermetic `aria-rnaseq-env` lockfile, `Dockerfile.rnaseq`, and C3 `bulk.preprocessing` ledger linkage |
| scRNA single-sample | Validated | PBMC 3k and GSE278576 single-sample runs |
| scRNA multi-sample | Validated | GSE278576 3-donor concat + Harmony workflow |
| Processed h5ad pseudobulk | Validated | 40-donor hippocampus rerun; report review confirmed pseudobulk, ORA, LIANA, trajectory, figures, and TSV exports |
| scATAC-seq universal raw & matrix engine (`scATAC`, ADR-076) | Validated (`production`) | Universal multi-platform raw FASTQ (`10x_atac`, `10x_multiome_atac`, `sci_atac`, `share_seq`, `surecell`, custom `--read-format`, multi-lane `L001..L00N`, automatic direct/RC barcode-whitelist orientation detection & fail-closed mismatch gate, `A4` resume), 100% offline air-gapped `.fai`/FASTA `snapatac2` fragment-to-matrix bridge, multi-sample fragments & `.h5mu`/`.h5ad` genomic-overlap consensus peak union (`ADR-043`), measured QC, TF-IDF/LSI, donor-aware `pyDESeq2` pseudobulk DA (`GSE278576` `5v5` `100%` subset precision vs R `DESeq2`), genomic peak annotation + local GMT ORA, `JASPAR2024` motifs, Signac-calibrated cis `peak2gene`, `ScATACAuditAgent` design governance, `Dockerfile.atacseq` + `Dockerfile.chromatin`, and `C3` `ledger://chromatin/scatac_preprocessing` linkage |
| Bulk ATAC-seq universal raw, BAM & matrix engine (`bulk_ATAC`, ADR-077) | Validated (`production`) | Universal SE/PE/multi-lane raw FASTQ (`fastp` + `atac_align.py` `bwa-mem2` with SE/PE flag geometry & `A4` resume), BAM/CRAM cohorts, and pre-computed peak-count TSV/CSV matrices; 100% local `pysam` + `bedtools` + local-GTF ENCODE QC (`FRiP`, `TSSe`, `dup_rate`, nucleosomal insert-size banding), SE/PE `MACS3` + consensus peaks (+ optional `blacklist_bed`), `bedtools` integer + `CPM` peak matrices (`chromatin_peak_counts.py`), replicate-aware `pyDESeq2` DA concordant with R `DESeq2`/`edgeR`/`limma-voom` (`docs/benchmark_results/bulk_atac/adr077_multireplicate_r_concordance.json`) + ENCODE `K562` vs `GM12878` motif recovery, genomic peak annotation + local GMT ORA, dual PNG+SVG PCA/QC/Volcano/MA figures (`chromatin_figure_bulk.py`), `BulkATACAuditAgent` replicate/confounding governance, and `C3` `ledger://chromatin/bulk_atac_preprocessing` linkage |
| Differential TF footprinting (`TOBIAS`, `scATAC` & `bulk_ATAC`, ADR-078) | Validated (`production_associative`) | Universal `TOBIAS 0.17.3` (`ATACorrect -> ScoreBigwig -> BINDetect` in locked `aria-tobias-env` / `Dockerfile.tobias`) across single- & multi-sample `scATAC` fragments and `bulk_ATAC` BAM/CRAM cohorts with automatic local FASTA/MEME resolution, `A4` `stage_manifest` resume, dual PNG+SVG `footprint_volcano` + aggregate profiles, live `footprint_rna_concordance`, `HypothesisAgent` `footprint_score` grounding (`FOOTPRINT_ASSOCIATIVE`), and external R concordance (`r = 1.0000`, `Jaccard = 1.000`, `0/20` null false discoveries across Welch, paired-donor, and batch-adjusted OLS in `docs/benchmark_results/bulk_atac_footprint/adr078_tobias_replicate_concordance.json`); pooled or `n<3` runs downgrade automatically to `descriptive_only` without FDR claims |
| Paired same-cell RNA+ATAC multiome engine (`INTEGRATION_VALIDATION`, ADR-079) | Validated (`production` / `production_associative`) | Promoted `INTEGRATION_VALIDATION` (`scope="paired_same_cell_rna_atac"`, `requires_ack=False`, `production_promoted=True`), all 8 integration `ScriptContract`s, `IntegrationAuditAgent`, and `Dockerfile.integration` (`ARIA_IMAGE_VALIDATION=validated`) to `production`: universal paired intake (`.h5mu`, Cell Ranger ARC `.h5` & dual MEX, shared/ARC-whitelist/bijection-mapped paired matrices, fragments+peaks, manifested raw Multiome FASTQs), sparse multi-donor consensus-peak cohort union (`integration_merge_cohort.py`), joint RNA+ATAC QC, donor/batch-aware Harmony + Seurat-concordant WNN (`ARI = 0.9582`) + 5-axis guardrail, ANOVA $\eta^2$ batch-screened MOFA+ (`|r| = 0.9730`), Signac-faithful cis `peak2gene` (`Jaccard = 0.8670`, `0` false discoveries on barcode-permuted multi-donor null control in `docs/benchmark_results/multiome_concordance/adr079_multidonor_multiome_validation.json`), candidate `eRegulons` (`17.99x` vs `SCENIC+`), `A4` `stage_manifest` resume across all stages, and `C3` `ledger://integration/preprocessing` linkage. Risk-elevated cohorts (`< 200` cells, partial barcode overlap, confounded donor/batch design) escalate to `yellow` `CP3.5` acknowledgement; unpaired/diagonal/bulk cross-modal integration stays blocked fail-closed |
| Integrated RNA biological discussion | Validated for RNA-only reports | Bulk RNA and scRNA synthesis blocks are generated from structured ARIA outputs, evidence-verified, claim-tiered, and rendered through the normal NarrativeAgent governance |

Current RNA/reporting baseline: release `v4.8.0` (building on `v4.7.0` with
scATAC completed to beta, bulk ATAC V47 beta lane, paired same-cell 10x Multiome
V48 beta integration + opt-in cross-modal `HypothesisAgent`, and the F1–F12 /
FASE 0–7 preprint-audit remediation — see below). The
post-tag hardening keeps donor-level pseudobulk as the primary scRNA inferential
layer, preserves raw QC counts for production pseudobulk, treats per-cluster
markers as descriptive rankings, removes prose-dependent QC thresholds, filters
marker candidates by `pct_in`, reuses compatible clustering marker rankings,
requires a defensible DPT root before pseudotime is computed, discloses
CellTypist fallback/low-confidence annotations, escalates destructive
integration-overcorrection signatures, treats a raw 10X MEX directory as one
sample, and renders failed QC as failed.

## Recent Hardening Closures

These changes explain why the `main` branch may be stricter than older reports:

| Closure | Status | What changed |
|---|---|---|
| scRNA W-CLAIM verification | Done | Diagnostic/error scRNA blocks, age-bin labels, thousands separators, LIANA tool names, and word-boundary analysis-family matching are handled without rejecting valid report text or accepting unsupported scientific claims |
| Anti-hardcode cleanup | Done | Runtime ligand-receptor fallback biology, MOFA cell-cycle mock biology, peak-to-gene mechanism labels from correlation sign, dataset/gene examples in production docstrings, and prose-derived CellTypist tissue hints were removed or gated |
| scRNA biological synthesis | Done | `BiologicalSynthesisAgent` now emits RNA-only scRNA integrated discussion blocks from measured pseudobulk, ORA, abundance, LIANA, trajectory, and reliability summaries; it remains data-only and makes no RNA+ATAC claims |
| Chromatin scATAC beta lane (ADR-048) | Done for beta | De-alpha'd in v4.7.0: scATAC QC/clustering, motif, and the replicate-gated pseudobulk condition-DA lane are beta-grade — externally concordant with SnapATAC2 (HC11 cluster ARI 0.532/NMI 0.669) and with edgeR-QLF/limma-voom/R-DESeq2 (5 young vs 5 old GSE278576 donors: ARIA ⊆ R-DESeq2 recall 1.000, LFC Spearman 0.61-0.76). Still dispatchable only behind explicit acknowledgement (`requires_ack`). The per-cluster Wilcoxon marker path stays caveated as single-sample-fragile; not publication-grade autonomy |
| Bulk ATAC replicate-gated DA | Beta opened | V47 comparison requests now build a peak-by-sample TSV from called/consensus peaks with `bedtools coverage -sorted -counts`, aggregate explicit biological replicates, and run the shared DESeq2 core over peaks when condition/replicate/comparison metadata are supplied. Missing metadata, absent contrasts, or insufficient replicates return structured skips; no filename-inferred contrasts or mock count matrices |
| Bulk ATAC TF motif interpretation | Beta opened | V47 DA peak sets are split by accessibility direction (both conditions, no one-sided pruning) and tested for hypergeometric TF motif over-representation against the tested-peak background, reusing the validated scATAC `chromatin_motifs` snapatac2 engine + versioned local JASPAR2024 collection (offline). Peaks on contigs absent from the reference FASTA are dropped before enrichment (disclosed). **Validated on real ENCODE K562 vs GM12878 DA peaks:** K562-up peaks recover the textbook erythroid KLF/SP signature (KLF1/EKLF, KLF5/7/15, SP1/2/4), GM12878-up peaks surface immune/B-cell IRF1 + Arid3a (`docs/benchmark_results/bulk_atac/v47_k562_gm12878_motif_enrichment.json`). Enrichment is associative (a database motif match), not evidence of TF activity/binding/regulation |
| scATAC beta extensions (P4) | Done for beta | On a real 10x PBMC Multiome: TOBIAS Tn5-bias-corrected TF footprinting + a **descriptive** BINDetect candidate ranking + RNA cross-evidence (single sample; associative, no causal or significance language), publication figures (dual PNG+SVG, TSS/FRiP gating panels), and gene-activity scoring (concordance vs Signac GeneActivity, moderate Spearman 0.51 → stays scaffold/caveated). B7 adds replicate-aware inference when an explicit barcode→group→replicate/donor design has ≥3 biological units per group: per-replicate mean-score tests, BH across motifs, paired-donor handling, and label-permutation diagnostics. Otherwise the lane stays descriptive. Dedicated `aria-tobias-env`; honest-skip when fragments/motifs/genome are absent. Tier remains beta/requires-ack. |
| Preprint-audit remediation (F1–F12) | Done | A 12-finding adversarial audit (Claude + Codex, cross-verified at file:line) fully remediated without touching the DE/DA numerical core: prompt-independent bulk-RNA `|log2FC|` cutoff (ADR-055); padj-ranked motif foreground; disclosed dropped DA covariates; no free-text LLM interpretation in Results; explicit bulk-RNA metadata + contrast contract; benchmark-artifact provenance; modality-correct motif thresholds; pydeseq2 fit warnings surfaced to the audit trail; disclosed public CSV-write failures; scoped V47 DA artifact + reproducible peak→gene marker mapping; machine-absolute-path scrub of committed graph artifacts; and a live synthetic recall/FDR gate for the bulk ATAC aggregation lane |
| Frozen open-weight reagent ablation (Track B3) | Done | Empirical 31-configuration hardware/determinism sweep across 5 open-weight families (`llama.cpp b12316` `5d4c8eeb2`, `Vulkan0`, `KV=f16`, `rb=0` → 25/25 configs `3/3` SHA-256 repeatable) + 18 live replays across ARIA's 6 core callsites (`CS1`–`CS6`) on pinned `Qwen3.6-35B-A3B-UD-Q4_K_M` (`100%` JSON validity, `0` ungrounded HGNC genes, `0` causal elevations, Claim 6 `Jaccard = 1.0000`) + `$0.00` repeated `run_headless` capsule equivalence ([supplement](benchmark_results/b3_frozen_reagent/B3_FROZEN_REAGENT_SUPPLEMENT.md)); hardened `aria/headless.py::_find_report_on_disk` so `ARIAMemory` never closes before post-HTML `memory_snapshot.sqlite` and `ro-crate-metadata.json` finish writing |
| Paired same-cell RNA+ATAC multiome integration & cross-modal hypotheses (V48 Slices I1–I9, ADR-064..ADR-072) | Beta opened (`requires_ack`) | `IntegrationAgent` promoted from `scaffold` to scoped **beta + `requires_ack`** (`scope="paired_same_cell_rna_atac"`, `production_promoted=False` per ADR-056): paired `.h5mu` joint RNA+ATAC QC (`I1`), WNN joint topology (`I2`) + 5-axis anti-degradation guardrail (`I3`), pure-Python MOFA+ variance decomposition with depth/cell-cycle confounder screening (`I4`), cis peak-to-gene GC/accessibility-matched background correlation (`I5`), candidate eRegulon triplet assembly with paralog-ambiguity & depth flags (`I6`), non-circular cross-modal evidence graph + `IntegrationNarrator` (`I7`), locked `aria-integration-env` (`I8`), and the opt-in `SPECULATIVE`-quarantined cross-modal `HypothesisAgent` adapter (`S10` / `I9`, `adapters/integration.py`, code-owned `independence_key` preventing cis double-counting). Unpaired, diagonal, and bulk cross-modal integration remain blocked fail-closed |
| Universal raw FASTQ RNA engine & `bulk_RNA_raw -> production` promotion (ADR-075) | Done (`production`) | Promoted `bulk_RNA_raw` and its 5 `rnaseq` IPC contracts (`rna_fastq_qc`, `rna_align`, `rna_umi_dedup`, `rna_quantify`, `rna_starsolo_bulk`) to `production`; added `rna_fastq_protocol.py` supporting `full_length`, `tag_3prime`, `tag_5prime`, `umi_bulk`, and `pooled_plate_bulk` (plus `scRNA` 3'/5' `kb` chemistry registry in `raw_ingestion.py`), fail-closed `pysam` multi-contig strandedness & A4 resume in `rna_quantify.py`, TPM length-bias suppression + CPM export for tag/UMI protocols, governed `bulk.preprocessing` C3 ledger linkage, `Dockerfile.rnaseq`, and locked `aria-rnaseq-env` (`pysam 0.23.3`, `numpy 1.26.4`) |
| Universal multi-platform scATAC engine & `scATAC -> production` promotion (ADR-076) | Done (`production`) | Promoted `scATAC` and its core spine IPC contracts (`chromatin_scatac_align`, `chromatin_fragments_to_matrix`, `chromatin_merge_matrices`, `chromatin_qc`, `chromatin_lsi_clustering`, `chromatin_diffacc`, `chromatin_peak_annotation`, `chromatin_peak_ora`, `chromatin_motifs`, `chromatin_figure_clusters`) to `production`; added multi-lane + multi-platform `chromap` geometries & auto-RC barcode orientation, offline air-gapped `.fai`/FASTA `snapatac2` matrix building, multi-sample fragments & `.h5mu`/`.h5ad` overlap consensus peak union, Signac-calibrated standalone `peak2gene`, `ScATACAuditAgent` design governance, `C3` `scatac_preprocessing` ledger linkage, and validated `Dockerfile.chromatin` + `Dockerfile.atacseq` |
| Universal bulk ATAC engine & `bulk_ATAC -> production` promotion (ADR-077) | Done (`production`) | Promoted `bulk_ATAC` (`BulkATACAuditAgent`) and its core IPC contracts (`atac_align`, `chromatin_peak_counts`, `chromatin_bulk_diffacc`, `chromatin_figure_bulk`, `chromatin_peak_annotation`, `chromatin_peak_ora`) to `production`; added SE/PE/multi-lane raw `bwa-mem2` alignment + `A4` resume, 100% local `pysam` + `bedtools` + local-GTF ENCODE QC (`FRiP`, `TSSe`, `dup_rate`, nucleosomal banding), SE/PE `MACS3` auto-detection + blacklist filtering, `CPM` peak-matrix export, pre-computed peak-count TSV intake, dual PNG+SVG PCA/QC figures, `C3` `bulk_atac_preprocessing` ledger linkage, and multi-replicate R `DESeq2`/`edgeR`/`limma-voom` concordance |
| Universal `TOBIAS` footprinting engine & `TOBIAS -> production` promotion (ADR-078) | Done (`production_associative`) | Promoted `chromatin_footprint_tobias.py`, `differential_tf_footprinting` (`scATAC` + `bulk_ATAC`), and `Dockerfile.tobias` (`ARIA_IMAGE_VALIDATION=validated`) to `production`; added multi-sample `fragments_files`, CRAM `--reference`, batch-adjusted OLS + stratified label permutations, `A4` `stage_manifest` resume, dual PNG+SVG `footprint_volcano`, auto-resolved FASTA/MEME + auto-derived `scATAC` peaks/design, live `footprint_rna_concordance`, `HypothesisAgent` `footprint_score` grounding (`FOOTPRINT_ASSOCIATIVE`), and external R + live `TOBIAS 0.17.3` concordance |
| Universal paired same-cell multiome engine & `INTEGRATION_VALIDATION -> production` promotion (ADR-079) | Done (`production`) | Promoted `INTEGRATION_VALIDATION` (`scope="paired_same_cell_rna_atac"`, `requires_ack=False`, `production_promoted=True`), all 8 integration IPC contracts (`integration_assemble_multiome`, `integration_merge_cohort`, `integration_paired_qc`, `integration_wnn`, `integration_wnn_diagnostics`, `integration_mofa`, `integration_peak2gene`, `integration_eregulons`), `IntegrationAuditAgent`, and `Dockerfile.integration` to `production`; added universal ARC `.h5`/MEX/paired-matrix/raw-FASTQ intake, multi-donor consensus-peak cohort union, donor/batch-aware Harmony WNN + WNN iLISI + overcorrection QC, MOFA+ ANOVA $\eta^2$ batch/donor screening, `A4` `stage_manifest` resume across all stages, `C3` `ledger://integration/preprocessing` provenance, and 4-donor `GSE278576` empirical + barcode-permuted null validation (`adr079_multidonor_multiome_validation.json`) |
| scRNA-lane production audit | Done | Closed B-PB1, B-DD1, B-QC1/B-QC2, A-MARK1, A-CMT1/A-CLUST1, and B-TRAJ1. Production pseudobulk uses raw QC counts, marker claims remain descriptive, QC is data-intrinsic, and DPT is skipped with `root_unresolved` when no defensible root exists |
| scRNA annotation/integration audit | Done | Closed N-ANNO1/N-ANNO2/N-ANNO3, N-QC1, and N-INT1. CellTypist confidence is genuine, low-confidence cell-type labels cap pseudobulk block confidence, default immune model fallback is disclosed, count-MAD QC is log-space, and destructive integration overcorrection is blocking |
| scRNA production E2E verification | Done | Real raw 10X pbmc3k runs verified the above report surfaces. Raw 10X MEX directories are collapsed to one sample, QC failures render as errors, and multi-donor pseudobulk verified `count_source=raw_counts` with the handoff |

## Reliability, Governance & Reproducibility

These cross-cutting guarantees back every validated path and are enforced in CI:

| Guarantee | Status | How |
|---|---|---|
| Deterministic narrative | Done | LLM calls at `temperature=0`; a fixed seed is passed only to backends that accept one, and each report declares whether the answering backend is seed-deterministic; model / tier / cache hit / degraded fallback recorded per report |
| Adversarial review of claims | Done | Deterministic, LLM-free devil's advocate enumerates the standard confounders (batch, ambient RNA, doublets, composition shift, low replication) per claim |
| Planned-vs-run ledger | Done | Report reconciles planned vs executed analyses; a partial run is visible, not silent |
| Single-source provenance | Done | version + git commit + dirty state + input SHA-256 + per-stage parameter hashes + dependency lockfiles + LLM usage |
| Design honesty | Done | confirmed covariates in the fitted DESeq2 formula; no alphabetical reference contrast; no filename-fallback design in production; user thresholds propagated end-to-end |
| No fabricated science | Done | repo-wide anti-fabrication guard (no placeholder matrices, ungated mock "successes", or hash-derived metrics) |
| No hardcoded biology in runtime fallbacks | Done | biological marker panels, ligand-receptor fallback tables, tissue keyword maps, and mechanism labels from feature names/signs are not used as production evidence |
| Evidence-linked claims | Done | report claims pass strict evidence verification, causal-language guards, run-ledger linkage, and deterministic devil's-advocate review before rendering |
| Untrusted prompt data + exclusive public claim path (C1) | Done | User questions, filenames, labels, genes, database terms and generated results enter LLM calls through one typed JSON `<untrusted_data>` boundary whose delimiters cannot be closed by payload text. Public HTML and `methodology.json` consume the same `compile_public_claims` output; unsupported blocks are withheld, and raw MessageBus/legacy prose is never promoted through an “All Findings” bypass |
| Semantic claim facts + explicit negation (C2) | Done | W-CLAIM carries serializable `SemanticFact` records with subject type, predicate, polarity, aliases, value and provenance. Differential-expression claims are checked by entity/alias and affirmation/negation, so a positive entity claim cannot be supported by “no significant genes”, a null claim cannot survive positive evidence, and Title-case symbols are resolved from the controlled predicate grammar rather than capitalization guesses. The prior numeric/entity checks remain defense-in-depth for predicates not yet modeled |
| Typed IPC contracts | Done | every dispatchable script validates inputs/outputs before and after the subprocess |
| Blocking CI (3 tiers) | Done | PR (guards + unit + contracts) / main (real pyDESeq2 recovery benchmark) / release (Docker env solve + in-image benchmark) |
| Air-gapped mode | Partial | Each execution now owns a live egress policy used by its LLM provider, pathway/connector gates, setup download guards, and child-process environment; two concurrent runs do not share that decision. `ARIA_AIR_GAPPED=1` remains the legacy launch-time default. Remaining gap (A1b): the initial biological-question LLM parse still occurs before CP1 sensitivity classification, so whole-run zero-egress cannot yet be claimed |
| Sensitivity checkpoint | Done | inputs are classified for clinical/PHI-like fields and quasi-identifiers before CP1; the user is always offered an air-gapped opt-in (never auto-disabled) |
| Local versioned enrichment | Done | over-representation analysis runs locally by default against versioned GMTs (library + release + SHA-256 recorded); Enrichr is opt-in and skipped honestly when egress is blocked |
| Deterministic packaging | Done | PEP 621 `pyproject.toml` with a single version source and per-platform lockfiles; `requirements.lock` core fallback |
| Hermetic build + supply chain | Done | per-modality Docker images with a `.dockerignore` (no `.git`/private `memory/`/caches in layers), gitleaks secret scan, and a CycloneDX SBOM of the RNA image; image digest stamped in the report |
| Secret hygiene | Done | installer reads API keys without echo; `aria doctor --secrets`/`--llm` classify/mask keys and flag committed credentials (no LLM call) |
| Per-modality readiness gate | Done | an assay capability matrix marks each modality green/yellow/red; green auto-dispatches, yellow needs explicit acknowledgement, red is removed from dispatch |

## Validated / Beta

| Area | Status | Notes |
|---|---|---|
| Raw-read ingestion (scRNA 3'/5' `kb count`) | Beta | Validated `kb` chemistry registry (`10XV1..10XV4`, `10X_5P` / `10X_5P_PE` TSO-aware geometry, `SMARTSEQ2/3`, `BDWTA`, `CELSEQ2`, `DROPSEQ`, `SPLIT-SEQ`, `STORMSEQ`, `VISIUM`, custom `bc:umi:cdna`) via `aria-ingestion-env` with SHA-256 reference and output provenance |
| Trajectory: PAGA + root-gated DPT | Beta | PAGA validated on hippocampus subset; DPT requires precomputed `iroot`, an explicit matching root label, or a generic progenitor/stem/precursor label. If no root is available, pseudotime is skipped with `root_unresolved`; trajectory remains exploratory, not causal |
| Cell-cell communication: LIANA | Beta | Validated on GSE278576 annotated h5ad; `n_perms=1000` default for stable ranks |
| GEO/SRA connector | Beta | **E6 atomic retrieval:** each accession is built in private staging, supplementary tar/tgz/zip members are traversal/link/size guarded, gzip/HDF5/MatrixMarket/FASTQ payloads are content-checked, every published byte is SHA-256 manifested, and an incomplete/corrupt generation never enters DataAudit or replaces a valid cache. When no analyzable GEO supplement exists, official NCBI RunInfo drives all-or-fail `prefetch` → `vdb-validate` → `fasterq-dump` per run; mixed-library FASTQs keep an authoritative per-file modality. Raw SRA retrieval requires the external SRA Toolkit, and downstream assay contracts still apply (for example scATAC needs its typed R1/R2/R3 + whitelist manifest). [Real-study evidence](benchmark_results/geo_connector/e6_atomic_retrieval.json) covers accession→DataAudit for all four modalities; live raw-SRA conversion remains outside that artifact. Public metadata remains heterogeneous; beta is not a promise that every accession contains a supported assay payload. |
| Raw-read ingestion (bulk ATAC) | Beta | Bulk ATAC FASTQ→BAM (`atac_align.py`, bwa-mem2) real-validated on ENCODE K562 chr20 (`bulk_ATAC` remains acknowledgement-gated `beta`; scATAC `chromap` + `snapatac2` raw ingestion is `production` under `ADR-076`). Honest-skip without the aligner / genome |
| scATAC optional exploratory sub-layers (`gene_scores` & per-cluster markers) | Caveated | While core `scATAC` ingestion, QC, TF-IDF/LSI, donor-aware `pyDESeq2` pseudobulk DA, peak annotation, local ORA, `JASPAR2024` motifs, Signac-calibrated `peak2gene`, and `TOBIAS` differential TF footprinting (`ADR-076`, `ADR-078`) are `production`, per-cluster Wilcoxon markers stay caveated as descriptive rankings, `gene_scores` stays caveated (`Spearman 0.51` vs Signac `GeneActivity`), and `TOBIAS` pooled/under-replicated (`n<3`) runs automatically downgrade to `descriptive_only` without FDR claims. |
| Paired multiome associative sub-layers (`peak2gene`, `eRegulons` & `MOFA+` factors) | Caveated (`production_associative`) | While paired same-cell RNA+ATAC intake, cohort union, joint QC, and Harmony/WNN topology are promoted to `production` (`ADR-079`), cis peak-to-gene correlations, candidate eRegulon triplets, and MOFA+ latent factors remain strictly observational (`production_associative`, zero causal verbs; WNN is reported as a joint topology lens alongside per-modality structure, never a replacement). |

## Alpha

No modalities are currently at alpha: scATAC was promoted to scoped **beta** in
v4.7.0 (ADR-048). Any future alpha modality dispatches only behind explicit
acknowledgement, remains review-required, and must surface missing resources, low
replication, and skipped lanes rather than fabricating output.

## Scaffolded / Roadmap

Scaffolded modalities are **dispatch-gated**: they cannot silently produce
publication-looking output. Direct calls to planned-but-absent scripts return a
structured `script_not_implemented`, and unpaired/diagonal cross-modal inputs
fail closed at the multiome contract gate rather than fabricating joint output.

| Area | Status | Required before stable |
|---|---|---|
| ChIP-seq / CUT&RUN / CUT&TAG | Scaffolded | Clear assay-specific QC and peak interpretation |
| Hi-C / Micro-C | Scaffolded — dispatch OFF | Runs only under `ARIA_ALLOW_EXPERIMENTAL_HIC=1` and are stamped not-publication-grade; needs E2E validation + memory-safe fixtures |
| Unpaired / diagonal / bulk cross-modal integration | Blocked fail-closed (deferred) | Confirmed same-cell RNA+ATAC integration is `production` (`ADR-079`, `scope="paired_same_cell_rna_atac"`); unpaired single-cell bridge/diagonal integration and bulk RNA+ATAC joint integration remain deferred |

Bulk ATAC-seq is no longer scaffolded: the V47 lane (QC + MACS3 peaks + peak×sample
counts + replicate-gated DESeq2 DA + TF-motif interpretation) is acknowledgement-gated
**beta**, validated end-to-end on real ENCODE replicates (see *Recent Hardening
Closures* above).

Chromatin QC (`chromatin_qc.py`) emits only measured metrics. TSS enrichment and
FRiP are real when their inputs/resources exist and otherwise remain null with a
concrete skip reason, never a fabricated placeholder.

The scATAC lane has been exercised on the local HC11 validation input
(`hc11_paired.h5mu`): ARIA reads ATAC modality `atac`, reports real dimensions
of 3,143 cells x 60,990 peaks, runs TF-IDF/LSI clustering, performs honest
single-sample DA where possible, and returns pseudobulk DA as skipped when
replicates/condition metadata are absent. The de-alpha to beta (v4.7.0 / ADR-048)
rests on external concordance: SnapATAC2 on HC11 (cluster ARI 0.532/NMI 0.669) and
a multi-replicate pseudobulk DA validation against edgeR-QLF/limma-voom/R-DESeq2 (5
young vs 5 old GSE278576 donors over a genomic-overlap consensus peak universe;
ARIA's 2,052 DA peaks ⊆ R-DESeq2's 6,646, LFC Spearman 0.61-0.76). The per-cluster
Wilcoxon marker path stays caveated as a descriptive ranking, while the core `scATAC`
modality (`chromap` universal raw FASTQ, air-gapped `snapatac2` fragments-to-matrix,
cohort consensus peak union, measured QC, TF-IDF/LSI, donor-aware `pyDESeq2` DA,
peak annotation, local ORA, `JASPAR2024` motifs, and Signac-calibrated `peak2gene`)
is promoted to `production` (`ADR-076`, governed by `ScATACAuditAgent`), and
`TOBIAS` differential TF footprinting is promoted to `production` (`ADR-078`,
`production_associative`), while `gene_scores` remains a caveated exploratory sub-layer. The workflow was subsequently
extended in P4 on a real 10x PBMC Multiome — peak-to-gene links (ADR-050),
TOBIAS footprinting + RNA cross-evidence, publication figures, and gene-activity
scoring — with gene-activity/motif-activity kept as caveated scaffold. B7 corrects
the footprint interpretation: that single-sample PBMC contrast is descriptive.
The replicate-aware path was separately exercised on real ENCODE K562 vs GM12878
biological replicates (2×2 explicit low-power override) over 7,374 chr21+chr22
peaks: 879 motifs tested, 209 BH calls, and zero discoveries in all four available
non-observed label permutations; the expected GM12878 IRF program (IRF8/4/9/5/2/1/7/3)
reappeared in the correct direction. This restricted-chromosome smoke validates
replicate plumbing, direction and the null-control contract, not genome-wide
sensitivity or a general FDR estimate; the deterministic synthetic null guard is
the calibration regression surface.

## Report Release Gate

Every generated report is held to these checks (originating in the `v4.3.12`
closeout and extended since):

- executive summary does not contradict body sections;
- pseudobulk DE, pathway, LIANA, and trajectory outputs are represented when present;
- missing outputs stay explicitly missing;
- supplementary TSV tables are non-empty when structured outputs exist;
- methods record the fitted design formula, thresholds, grouping columns, covariates, and warnings;
- existing `obs` annotations are reported as reused groupings, not as newly inferred Leiden clusters;
- claims are evidence-tiered and never exceed the language their evidence licenses;
- every rendered scientific claim must be backed by local evidence cards or be
  withheld; diagnostic blocks can report structured skip/error counts without
  masquerading as measured biology;
- modeled claim predicates must match typed semantic evidence in entity/alias and
  polarity; explicit null evidence cannot license an affirmative result sentence;
- the submitted question is visibly labeled untrusted input and is not evidence
  for an executive-summary claim; raw MessageBus text remains an operational audit
  record, while the public claim ledger is generated only from compiled
  `NarrativeBlock` manifests;
- `rna_de_per_cluster.py` is treated as optional on atlas-scale inputs. If it times
  out, the pipeline may still be scientifically valid when donor-level pseudobulk,
  pathway, communication, and trajectory outputs complete. DPT pseudotime may be
  absent by design when no defensible root is available.

Latest reviewed hippocampus rerun produced real output tables and figures; newer
reports state that Leiden was skipped when input annotations are reused, and carry
the full provenance and evidence-tier manifest described above. The integrated
RNA discussion is part of the governed report path; paired same-cell RNA+ATAC
multiome synthesis is available in the `production` paired-multiome lane (`IntegrationNarrator`,
ADR-070/ADR-079, risk-gated by `IntegrationAuditAgent`), while unpaired cross-modal synthesis remains blocked fail-closed.
