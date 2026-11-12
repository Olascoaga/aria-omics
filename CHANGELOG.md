# Changelog

All notable milestone releases of **ARIA (`aria-omics`)** are documented here.
`aria/version.py` is the single runtime version source; historical release notes
live under [`docs/releases/`](docs/releases/) and describe the state at each past
tag, whereas [`docs/validation_status.md`](docs/validation_status.md) reflects
current `main`.

---

## [Unreleased] — `main` (post-`v4.8.0`)

### Highlights
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
