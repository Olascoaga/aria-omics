# Changelog

All notable milestone releases of **ARIA (`aria-omics`)** are documented here.
`aria/version.py` is the single runtime version source; historical release notes
live under [`docs/releases/`](docs/releases/) and describe the state at each past
tag, whereas [`docs/validation_status.md`](docs/validation_status.md) reflects
current `main`.

---

## [Unreleased] — `main` (post-`v4.8.0`)

### Highlights
- **Universal `TOBIAS` Differential TF Footprinting Engine & `TOBIAS -> production` Promotion (`ADR-078`):**
  - Promoted `aria/scripts/chromatin_footprint_tobias.py` (`ScriptContract`), `differential_tf_footprinting` (`production_associative` across `ScATACAuditAgent` and `BulkATACAuditAgent`), and `docker/Dockerfile.tobias` (`ARIA_IMAGE_VALIDATION=validated`) to **`production`** (`INTEGRATION_VALIDATION` remains scoped `beta + requires_ack`).
  - Upgraded `chromatin_footprint_tobias.py` with multi-sample `fragments_files` streaming (`sample#barcode` / `_` / `:` prefix resolution), CRAM `--reference` support, additive batch-adjusted linear modeling (`score ~ 1 + batch + group` + stratified within-batch label permutations + `batch_condition_confounded` honest downgrade), cryptographic `A4` `stage_manifest` resume, and dual PNG+SVG `footprint_volcano` rendering.
  - Wired `ChromatinAgent` for zero-friction local `genome_fasta`/`.fai` and `JASPAR2024` MEME motif auto-resolution, auto-derived `scATAC` `peaks.bed` & `barcode_design.tsv` from clustered `.h5ad`, raw-FASTQ `sample_ids` & `replicate_batches` propagation in `bulk_ATAC`, live `footprint_rna_concordance` attachment, copy-pasteable `ChromatinNarrator` Methods + Tables, `HypothesisAgent` `footprint_score` grounding (`FOOTPRINT_ASSOCIATIVE`, refusing `descriptive_only` runs), and external R (`t.test` Welch & paired + `lm` batch OLS + BH) + live `TOBIAS 0.17.3` concordance (`docs/benchmark_results/bulk_atac_footprint/adr078_tobias_replicate_concordance.json`).
- **Universal `Bulk ATAC-seq` Production Engine & `bulk_ATAC -> production` Promotion (`ADR-077`):**
  - Promoted `MODALITY_VALIDATION["bulk_ATAC"]` (`BulkATACAuditAgent`: `green`/`auto` for non-contrasted profiling or `n >= 3` biological-replicate contrasts, `yellow`/`requires_ack` for `n == 2` low-power isogenic replicates, `red`/`blocked` for `n < 2`, batch-confounded contrasts, or single-cell barcoded FASTQs misrouted to `bulk_ATAC`) and the `bulk_ATAC` IPC spine (`atac_align.py`, `chromatin_peak_counts.py`, `chromatin_bulk_diffacc.py`, `chromatin_figure_bulk.py`, `chromatin_peak_annotation.py`, `chromatin_peak_ora.py`) to **`production`** (`TOBIAS` footprinting and `INTEGRATION_VALIDATION` remain scoped `beta + requires_ack`).
  - Upgraded `atac_align.py` with full Paired-End and Single-End (`paired=False` omitting `fixmate`/`-f 2` while retaining `markdup -r`, MAPQ, and mito stripping) + multi-lane / `bulk_fastq_manifest` intake and `A4` stage-manifest resume; hardened `chromatin_qc.py` and `chromatin_peaks.py` with 100% local `pysam` + `bedtools` + local-GTF ENCODE QC (`FRiP`, `TSSe`, `dup_rate`, nucleosomal banding without `samtools` in `aria-chromatin-env`), SE/PE BAM auto-detection (`--nomodel --shift -100 --extsize 200`), optional `blacklist_bed` filtering, and removal of the 8-sample QC cap.
  - Promoted `chromatin_peak_counts.py` (`CPM` matrix export + CRAM + `A4` resume), fixed `_run_bulk_atac` contrast gating so `design.contrasts` / `comparisons` trigger `peak_counts -> DA`, added direct pre-computed peak-count TSV/CSV intake, added `chromatin_figure_bulk.py` (`bulk_atac_pca` + `bulk_atac_qc_summary` dual PNG+SVG), closed `C3` `ledger://chromatin/bulk_atac_preprocessing` linkage, and verified multi-replicate (`n=4v4` + technical replicates + batch covariate + 20 null permutations) external concordance against R `DESeq2`, `edgeR`, and `limma-voom`.
- **Universal Multi-Platform `scATAC-seq` Engine & `scATAC -> production` Promotion (`ADR-076`):**
  - Promoted `MODALITY_VALIDATION["scATAC"]` (`ScATACAuditAgent`) and the core `scATAC` IPC spine (`chromatin_scatac_align.py`, `chromatin_fragments_to_matrix.py`, `chromatin_merge_matrices.py`, `chromatin_qc.py`, `chromatin_peaks.py`, `chromatin_lsi_clustering.py`, `chromatin_diffacc.py`, `chromatin_peak_annotation.py`, `chromatin_peak_ora.py`, `chromatin_motifs.py`, `chromatin_figure_clusters.py`) to **`production`**.
  - Upgraded `scatac_fastq_manifest.py` and `chromatin_scatac_align.py` with multi-lane (`L001..L00N`) FASTQ lists per role, platform presets (`10x_atac`, `10x_multiome_atac`, `sci_atac`, `share_seq`, `surecell`, custom `--read-format`), automatic direct vs reverse-complement barcode-whitelist orientation detection (failing closed on `barcode_whitelist_mismatch`), and `A4` `stage_manifest` resume.
  - Hardened `chromatin_fragments_to_matrix.py` for 100% offline `ARIA_AIR_GAPPED=1` matrix building from local `.fai`/FASTA/`chrom_sizes`, multi-sample `fragments_files` union, `obs` metadata preservation, and `A4` resume; added `scatac_cohort_union.py` + `chromatin_merge_matrices.py` for sweep-line genomic-overlap consensus peak union across multi-file `.h5mu`/`.h5ad` cohorts (`ADR-043`).
  - Closed `scATAC` analytical parity with per-comparison `full_results_csv` in `chromatin_diffacc.py`, `chromatin_peak_annotation.py` + `chromatin_peak_ora.py` wiring in `ChromatinAgent`, Signac-calibrated `(GC, accessibility)` trans-null `peak2gene` in `chromatin_regulatory.py`, `C3` `ledger://chromatin/scatac_preprocessing` narrative linkage, and validated `docker/Dockerfile.chromatin` + `docker/Dockerfile.atacseq` + `beta` `docker/Dockerfile.tobias`.
- **Universal Raw FASTQ Multi-Platform RNA Engine & `bulk_RNA_raw -> production` Promotion (`ADR-075`):**
  - Promoted `MODALITY_VALIDATION["bulk_RNA_raw"]` and all five `aria-rnaseq-env` IPC contracts (`rna_fastq_qc.py`, `rna_align.py`, `rna_umi_dedup.py`, `rna_quantify.py`, `rna_starsolo_bulk.py`) to **`production`**.
  - Added `aria/utils/rna_fastq_protocol.py` supporting four bulk RNA-seq platform families from raw FASTQ: `full_length` (SE/PE/multi-lane `STAR` 2-pass), `tag_3prime` & `tag_5prime` (`QuantSeq`, `SMART-3Seq`, `CAGE`/`RAMPAGE` with poly(X)/poly(G) trimming and 1-pass soft-clip-preserving alignment), `umi_bulk` (`fastp --umi` + `rna_umi_dedup.py` directional 1-Hamming `pysam` UMI graph collapse), and `pooled_plate_bulk` (`BRB-seq`, `DRUG-seq`, `Prime-seq`, `QuantSeq-Pool`, `CEL-Seq2` via `rna_starsolo_bulk.py` `STAR --soloType CB_UMI_Simple --soloUMIdedup 1MM_Directional_UMItools`).
  - Expanded `scRNA` raw FASTQ `kb-python` chemistry resolution (`aria/utils/raw_ingestion.py`) across 3', 5' (`10X_5P`, `10X_5P_PE`), full-length/UMI plate (`SMARTSEQ2`, `SMARTSEQ3`, `CELSEQ2`), combinational (`SPLIT-SEQ`, `BDWTA`), and custom `bc:umi:cdna` geometries, while failing closed on long-read (`>=600 bp`), V(D)J TCR/BCR repertoire, and unmanifested asymmetric barcode+UMI FASTQs.
  - Hardened `rna_quantify.py` with `A4` `counts_matrix.stage.json` resume, indexed/unindexed multi-contig `pysam` fail-closed strandedness detection (`AmbiguousStrandedness` / `InsufficientReadsForStrandDetection`), `CPM` export (`counts_cpm.tsv`) with automatic length-bias `TPM` suppression on tag/UMI protocols, `C3` `ledger://bulk/preprocessing` narrative linkage, locked `pysam>=0.22` in `aria-rnaseq-env`, and `docker/Dockerfile.rnaseq`.

---

## [v4.8.0] — 2026-11-12 (Current Release)

Full release notes: [`docs/release_notes_v4.8.0.md`](docs/release_notes_v4.8.0.md)

### Highlights
- **Paired Same-Cell 10x Multiome Integration (`beta + requires_ack`, V48 Slices `I1`–`I8`, `ADR-064`..`ADR-071`):**
  - Fail-closed `enforce_same_cell_integration_gate`, `read_h5mu_paired`, `h5py` `/obs`-only bijective label join, and joint RNA+ATAC QC (`integration_paired_qc.py`).
  - Deterministic Weighted Nearest Neighbors (`integration_wnn.py`, Seurat v5 `ARI = 0.9582`, `NMI = 0.9660`) paired with a 5-axis anti-degradation guardrail (`wnn_diagnostics.py`, "lens, not replacement").
  - Deterministic paired MOFA+ (`integration_mofa.py`, `muon.tl.mofa` `R^2`-weighted `|r| = 0.9730`) with closed-form variational subsample projection and automatic technical depth/cell-cycle factor exclusion.
  - Unified Signac-faithful cis `peak2gene` engine (`peak2gene_core.py` + `integration_peak2gene.py`, `Signac::LinkPeaks` `Jaccard = 0.8670`, `score_pearson = 1.0000`).
  - Deterministic paired candidate eRegulons (`eregulon_core.py` + `integration_eregulons.py`, `17.99x` hypergeometric edge fold-enrichment vs `SCENIC+ 1.0.0`).
  - Non-circular Cross-Modal Evidence Graph (`cross_modal_triangulation.py`), `IntegrationNarrator`, 4 dual PNG (`300 DPI`) + editable SVG publication figures, CP3.5 acknowledgement gate, and hermetic `aria-integration-env` Conda + PyPI lockfiles.
- **Cross-Modal Opt-In `HypothesisAgent` (`S10`, V48 Slice `I9`, `ADR-072`):**
  - Added 5th deterministic evidence adapter (`adapters/integration.py`), code-owned `independence_key` non-circular cis lineage ranking, and 4 scoped multiome caveats behind the unmodified 4-gate `SPECULATIVE` quarantine wall (`assert_no_speculative_promotion`).
- **Preprint Governance & Reproducibility Hardening:**
  - Closed `C1`–`C4` public-claim compilation, `SemanticFact` polarity verification, exact per-estimand causal licensing, `E1`–`E6` raw/GEO intake boundaries, `A2`/`A3`/`A5`/`A6` execution isolation and portable `capsule.zip` / RO-Crate exports, Track `B3` frozen open-weight reagent ablation, and the 21-lane `preprint_v1` cryptographic evidence freeze inventory (`18/21` machine lanes verified).

---

## [v4.7.0] — 2026-07 (`docs/releases/release_notes_v4.7.0.md`)
- Textual **ARIA Control Center** (`aria/ui/cockpit.py`, `intake.py`, `render.py`) over a pure execution read-model.
- End-to-end **Bulk ATAC-seq `beta + requires_ack`** workflow (`bwa-mem2` alignment, measured QC, `MACS3` consensus peaks, `bedtools` peak counts, replicate-gated `pyDESeq2` DA, genomic annotation, ORA, motif enrichment, and `TOBIAS` footprinting).
- Single-modality opt-in **`HypothesisAgent` (`S0`–`S9`, `ADR-057`)** with 4 deterministic publication gates and non-promotable `SPECULATIVE` quarantine.
- De-alpha promotion of **scATAC-seq to `beta + requires_ack` (`ADR-048`)** backed by multi-donor `GSE278576` pseudobulk DA concordance against R `DESeq2`, `edgeR-QLF`, and `limma-voom` plus external `SnapATAC2` and `Signac` concordance suites.

---

## [v4.6.x] — 2026-06 (`docs/releases/release_notes_v4.6*.md`)
- **[v4.6.1](docs/releases/release_notes_v4.6.1.md) & [v4.6](docs/releases/release_notes_v4.6.md):** scATAC raw FASTQ (`chromap`) and fragments-to-peak-matrix (`snapatac2`) bridge, `TOBIAS` footprinting integration, `Signac` peak-to-gene / gene-activity concordance harnesses, and multi-sample consensus peak unification (`ADR-043`).
- **[v4.6.0-alpha](docs/releases/release_notes_v4.6.0-alpha.md):** Initial live-orchestrated scATAC matrix pipeline (`QC`, `TF-IDF/LSI`, `diffacc`, `JASPAR2024` motif enrichment) behind CP3.5 readiness acknowledgement (`ADR-033`, `ADR-034`).

---

## [v4.5.x] — 2026-06 (`docs/releases/release_notes_v4.5.*.md`)
- **[v4.5.5](docs/releases/release_notes_v4.5.5.md):** Frozen v4.5 RNA & governance benchmarking spine (`A1` synthetic + `SEQC/MAQC` TaqMan + 5-site + `ERCC`, `A2` synthetic donor-aware + `Kang/muscat`, `B1` adversarial design refusal, `B2` 4-arm claim faithfulness, `B4` null-control narrative).
- **[v4.5.4](docs/releases/release_notes_v4.5.4.md) & [v4.5.2](docs/releases/release_notes_v4.5.2.md):** Guided scRNA raw FASTQ ingestion (`kb-python`), `BiologicalSynthesisAgent` (`ADR-028`), `W-LEDGER` claim-to-ledger verification (`ADR-026`), PEP 621 packaging, Conda lockfiles, and per-modality Dockerfiles.

---

## [v4.3.12 – v4.3.19] — Historical Archive
- Early RNA orchestration, design checkpoints, parameter advisor, and narrative hardening notes archived in [`docs/releases/`](docs/releases/).
