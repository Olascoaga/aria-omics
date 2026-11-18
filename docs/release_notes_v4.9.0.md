# ARIA v4.9.0 Release Notes

`v4.9.0` is a **minor release** over `v4.8.0`. It completes the promotion of **all
`9 / 9` single-assay modalities and all `3` cross-cutting multi-modal engines to
`production`** (`ADR-075`..`ADR-082`), introduces **pre-execution cryptographic
`DecisionManifest`s, append-only `AnalyticalForkLedger`s, Four-Layer Determinism,
and Zero-LLM Replay** (`ADR-084`), closes **walk-away raw `FASTQ` reference
auto-provisioning parity** across single-cell and paired Multiome workflows (`Q1`),
wires **live `In-Silico Reviewer #2` biological-unit LOBO jackknife falsification
certificates** across all four production DE/DA engines (`P1`), and ships the
**`ADR-083` Automated Preprint `v1` Gate (`19 / 19` machine lanes)**, the **600-case
metamorphic report falsification suite**, the **40-capsule governed vs.
uncontracted `CodeAct` comparative benchmark**, the **11-Modality `$0.00` Public
Showcase Atlas**, and the **30-Cell Cross-Model Equivalence & Cost Benchmark**.
The **`pyDESeq2` statistical core (`aria/scripts/rna_bulk/deseq2.py`) and
`aria/utils/stage_manifest.py` have zero diff**.

## Headline 1: Universal Production Analytical Spine (`9 / 9` Modalities + `3` Multi-Modal Engines) — ADR-075..082

All nine single-assay modalities in `orchestrator_agent.MODALITY_VALIDATION`
(`bulk_RNA`, `bulk_RNA_raw`, `scRNA`, `scATAC`, `bulk_ATAC`, `ChIP`, `CUT_AND_RUN`,
`CUT_AND_TAG`, `HiC`) and all three cross-cutting/multi-modal engines (`TOBIAS`
differential TF footprinting `ADR-078`, `INTEGRATION_VALIDATION` paired same-cell
Multiome `ADR-079`, and `CROSS_MODAL_CONCORDANCE_VALIDATION` non-same-cell RNA
$\times$ ATAC contrast & donor concordance `ADR-080`) now operate at
**`production`** across **9/9 lock-audited Conda stacks** (`environment_audit == "pass"`):

- **Universal Raw FASTQ Multi-Platform RNA (`ADR-075`):** `bulk_RNA_raw` supports
  `full_length` (`STAR` 2-pass), `tag_3prime`/`tag_5prime` (`QuantSeq`, `SMART-3Seq`,
  `CAGE`), `umi_bulk` (`fastp --umi` + directional 1-MM `pysam` dedup), and
  `pooled_plate_bulk` (`BRB-seq`, `DRUG-seq`, `Prime-seq`, `CEL-Seq2` via
  `STARsolo`), plus expanded `kb-python` `scRNA` chemistries (`10x` 3'/5',
  `Smart-seq2/3`, `BD Rhapsody`, `CEL-Seq2`, `SPLiT-seq`, custom `bc:umi:cdna`).
- **Universal Multi-Platform `scATAC`, `bulk_ATAC` & `TOBIAS` (`ADR-076`..`078`):**
  multi-lane `chromap` presets (`10x_atac`, `10x_multiome_atac`, `sci_atac`,
  `share_seq`, `surecell`, auto-RC barcode orientation), offline air-gapped
  `snapatac2` matrix construction, sweep-line consensus peak union, SE/PE `MACS3` +
  `bedtools` integer/CPM matrices, replicate-gated `pyDESeq2` DA, and inferential
  multi-sample `TOBIAS 0.17.3` differential TF footprinting (`production_associative`).
- **Universal Epigenomics (`ChIP`, `CUT_AND_RUN`, `CUT_AND_TAG`, `ADR-081`):**
  assay-aware `bwa-mem2` alignment (`remove_duplicates=False` default for
  `CUT_AND_TAG`) + heterologous spike-in (`sacCer3` / `E. coli`) read counting,
  ENCODE `NRF`/`PBC1`/`PBC2`/`NSC`/`RSC`/`JSD` QC, Meers et al. (2019) `SEACR`
  (`stringent`/`relaxed` IgG & `top-q`), broad megadomains, Li et al. (2011)
  Gaussian-copula EM `IDR`, BigWig export, and technical-replicate sum-first
  spike-in `size_factors` `pyDESeq2` differential binding.
- **Universal 3D Genome Architecture (`HiC` / `Micro-C`, `ADR-082`):**
  `GenomeArchAgent` + locked `aria-hic-env` (`cooler 0.10.4`, `cooltools 0.7.1`,
  `pairtools 1.1.3`, `hic2cool 1.0.1`) across `.mcool`, `.cool`, `.hic`,
  `.pairs(.gz)`, chimeric BAM, and raw paired FASTQ (`bwa-mem2 -SP5M`), with ICE
  balancing, $P(s)$ decay, HiCRep `SCC`, GC-phased `E1` compartments, Crane
  insulation TADs, `dotfinder` loops, and empirical-Bayes moderated `differential_3d`.
- **Paired Same-Cell Multiome (`ADR-079`) & Non-Same-Cell Cross-Modal Concordance
  (`ADR-080`):** universal paired `.h5mu` / Cell Ranger ARC `.h5`/MEX / raw Multiome
  FASTQ intake + Harmony-adjusted WNN + batch-screened MOFA+ + `peak2gene` +
  `eRegulons`, complemented by `cross_modal_concordance.py` for aligned `bulk_RNA +
  bulk_ATAC` and unpaired cohort `scRNA + scATAC` studies (while keeping diagonal
  single-cell manifold imputation fail-closed).

## Headline 2: Pre-Execution `DecisionManifest`, `AnalyticalForkLedger`, Four-Layer Determinism & Zero-LLM Replay — ADR-084

- **Committed `DecisionManifest` & Forking-Path Ledger (`aria/utils/decision_manifest.py`):**
  before any modality compute stage runs, `OrchestratorAgent` signs and writes
  `decision_manifest.json` (SHA-256 over canonical study input digests, design
  contract, contrasts, FDR/LFC thresholds, and parameter policy) and appends an
  immutable record to per-study `fork_ledger.jsonl` (preserving `fork_index=0` on
  exact replays while stamping `post_hoc_fork` and `fork_changes` whenever
  thresholds, covariates, or contrasts diverge on the same `study_input_fingerprint`).
- **Zero-LLM Replay (`run_headless(..., replay_manifest=...)`):** replays any
  verified `decision_manifest.json` end-to-end with **`0` LLM calls** across
  `OrchestratorAgent`, `DesignAgent`, `ParameterAdvisor`, and `NarrativeAgent`, and
  verifies `decision_manifest.json` inside `capsule.zip`.
- **Four-Layer Determinism & Tamper-Checked Stage Caches (`A4` / `A6`):**
  session-scoped BLAS/OpenMP/Numba single-thread pinning under `reproducible_mode` +
  `PYTHONHASHSEED=0` (`privacy.py`), `bwa-mem2 mem -K 100000000` (`atac_align.py`),
  explicit LIANA seed + `mergesort` tie-breaking (`rna_cellcomm.py`), deterministic
  ARPACK `v0` initialization (`hic_topology.py`), token-boundary keyword matching
  (`text_matching.py`), canonical 12-point Leiden grid `[0.20..1.40]`
  (`parameter_advisor.py`), and output SHA-256 tamper verification + advisory file
  locks in `mudata_io._stage_cache`, `hic_preprocess`, and `rna_advise_resolution`.

## Headline 3: Walk-Away Raw `FASTQ` Reference Auto-Provisioning Parity (`Q1`) — ADR-084

- **Lock-Serialized `kb ref` & `chromap.index` Provisioning:** `SetupAgent._ensure_kb_index`
  automatically builds and SHA-256-verifies `kb-python` transcriptome indices
  (`index.idx` + `t2g.txt`) from staged `FASTA + GTF` under advisory `fcntl` locks,
  consumed automatically by `RawIngestionAgent` and `IntegrationAgent`;
  `integration_assemble_multiome.py` likewise auto-builds `chromap.index` when
  absent, and `propose_illumina_scatac_manifest` discovers standard Illumina 10x
  ATAC folder layouts (`R1`/`R2`/`R3`/`I1`) for CP1 confirmation. Sequencer-native
  pre-demultiplex/basecalling formats (`BCL`, `POD5`, `FAST5`) remain fail-closed (`E1`).

## Headline 4: Live `In-Silico Reviewer #2` Biological-Unit Falsification Certificates (`P1`)

- **Stack-Isolated LOBO Engine (`aria/utils/lobo_jackknife.py`):** implements
  Smyth (2004) empirical-Bayes variance-moderated Leave-One-Biological-Unit-Out
  jackknife across both unpaired biological replicates and paired-donor designs
  (preserving DESeq2 median-of-ratios or heterologous spike-in `size_factors` and
  refusing duplicate technical-replicate unit IDs or `n < 3` arms).
- **Live DE/DA Integration & Orthogonal Certificates:** `rna_bulk_de.py`,
  `rna_pseudobulk_de.py`, `chromatin_diffacc.py`, and `chromatin_bulk_diffacc.py`
  attach live contrast-level `lobo_jackknife` diagnostics on every `n >= 3` fit,
  automatically extracted into 4-axis `FalsificationCertificate`s
  (`falsification_certificates.py`) in `methodology.json` and `report.html`
  (`feeds_into_parameter_selection=False`, zero scalar fragility score, zero causal
  promotion).

## Headline 5: `ADR-083` Automated Preprint Gate, 600-Case Metamorphic Suite, 40-Capsule Comparative Benchmark & Public Vignettes `01`–`13`

- **Automated Preprint `v1` Gate & 600-Case Metamorphic Suite (`ADR-083`):**
  `preprint_freeze.py` exposes `automated_v1_gate` (`19 / 19` machine-verifiable
  lanes across `Claims 1–7`, keeping the `3` `human_gold` reader-study lanes
  untouched at `0` diff as `deferred_to_v2_reader_study`), including
  `c4_metamorphic_600_falsification` (`300 / 300` invalid adversarial mutations
  caught across 6 error families $\times$ 5 modalities and `300 / 300` valid
  controls retained, paired with tightened multi-omics negation polarity in
  `evidence_verifier.py`, `failed`/`blocked`/`cancelled` status rejection in
  `run_ledger.py`, and `"cause"`/`"causally"` lexicon closure in `validators.py`).
- **40-Capsule Governed vs. Uncontracted `CodeAct` Benchmark (`codeact_comparative_suite.py`,
  Vignette `13`):** evaluates `aria_governed`, `aria_locked_manifest`,
  `aria_zero_llm_deterministic`, and `uncontracted_codeact_baseline` via
  `aria/benchmarks/adapter.py` + `derived_query_ast.py` across 40 deterministic
  capsules (`5` archetypes $\times$ `8` capsules): `aria_governed` achieves
  `silent_error_rate = 0.000` (`1.000` conditional accuracy, `1.000` abstention P/R,
  `40/40` locked-manifest & zero-LLM concordance; `0.0025` mean silent-error rate
  and `0` trap violations across a 20-seed `800`-run sweep) vs `0.650` (`26/40`) on
  the uncontracted baseline.
- **11-Modality `$0.00` Public Showcase Atlas (Vignettes `01`–`11`) & 30-Cell
  Cross-Model Equivalence Benchmark (Vignette `12`):** full peer-reviewed public
  GEO/ENCODE/ENA capsules and multi-provider cost/invariance receipts committed
  under `docs/vignettes/` and `docs/benchmark_results/`.
