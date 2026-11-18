# ARIA v4.8.0 Release Notes

`v4.8.0` is a **minor release** over `v4.7.0`. It delivers the **paired same-cell
10x Multiome (RNA + ATAC) integration stack** (`ADR-064`..`ADR-071`, promoted from
`scaffold` to scoped `beta` + `requires_ack`), the opt-in **`HypothesisAgent`
(`S0`–`S10`)** including its 5th cross-modal evidence adapter (`ADR-057`..`ADR-062`,
`ADR-072`), the **bulk ATAC-seq `beta` lane** (`ADR-054`), the **F1–F12 / FASE 0–7
preprint-readiness governance hardening** (`C1`–`C4`, `A1`–`A7`, `B1`–`B7`,
`E1`–`E7`), and the **Track B3 frozen open-weight reagent ablation**. The **DE/DA
estimation core is unchanged**. This is a normal package release, **not** the
immutable `preprint-v1` artifact freeze tag (the three independent human-gold
design/claim adjudication lanes remain explicit freeze blockers).

## Headline 1: Paired Same-Cell 10x Multiome Integration (`beta` + `requires_ack`) — ADR-064..ADR-071

`orchestrator_agent.INTEGRATION_VALIDATION` is promoted from `scaffold` to scoped
**`beta` + `requires_ack=True`** (`scope="paired_same_cell_rna_atac"`,
`production_promoted=False` per `ADR-056`), while `MODALITY_VALIDATION` stays
strictly at its 9 single-assay keys and `aria-integration-env` is promoted to
`lock_required=True` (`envs/aria-integration-env.{linux-64.lock,pip.lock}`).

- **Fail-closed scope & CP3.5 gate (`I1`, `I8`):** `enforce_same_cell_integration_gate`
  requires confirmed shared-barcode same-cell pairing (`.h5mu` or aligned same-cell
  objects) across every `IntegrationAgent` strategy; unpaired split `scRNA + scATAC`,
  multi-`.h5mu` directories without an explicit merge map, and bulk cross-assay
  combinations fail closed without blocking independent single-assay lanes. Dispatch
  additionally requires an affirmative CP3.5 acknowledgement logged in `ARIAMemory`.
- **Joint paired QC & bijective reference transfer (`I1`, ADR-064):**
  `integration_paired_qc.py` filters shared barcodes across raw RNA + ATAC counts
  and supports `h5py` `/obs`-only bijective barcode joins (`2,517/3,143` post-QC
  cells on `GSE278576` `hc11_paired.h5mu`, `100.0%` 1-to-1 `subclass` match).
- **WNN joint topology + 5-axis anti-degradation guardrail (`I2`–`I3`, ADR-065..066):**
  deterministic `integration_wnn.py` externally concordant with R `Seurat 5.3.0`
  `FindMultiModalNeighbors` on `hc11` (`ARI = 0.9582`, `NMI = 0.9660`, modality-weight
  `Spearman rho = 0.9277`, `MAE = 0.0275`), paired with `wnn_diagnostics.py` verifying
  per-modality recall, rare-subpopulation retention (`6/6` evaluable rare subclasses
  retained), zero spurious bridges, and out-of-fold cross-modal cosine lift (`0.5786`
  joint vs `0.5204` RNA-only and `0.4293` ATAC-only). WNN is governed strictly as an
  additional lens, never a replacement for per-modality structure.
- **Deterministic MOFA+ variance decomposition (`I4`, ADR-067):** `integration_mofa.py`
  with Horvitz-Thompson stratified subsampling, closed-form Gaussian variational
  posterior projection (`0 NaN`s), hybrid sparse-view activity classification, and
  automatic sequencing-depth (`|r| >= 0.70`) / cell-cycle confounder exclusion;
  externally concordant with `muon.tl.mofa` on `hc11` (Hungarian `R^2`-weighted
  `|r| = 0.9730`, canonical subspace correlation `0.9895`).
- **Unified cis peak-to-gene engine (`I5`, ADR-068):** shared `peak2gene_core.py`
  with bidirectional UCSC/Ensembl contig canonicalization, `.fai` byte-offset GC
  computation, and vectorized `(GC, accessibility)` trans-chromosomal background null;
  externally concordant with `Signac::LinkPeaks` on `hc11` (`Jaccard = 0.8670`,
  `score_pearson = 1.0000`, `zscore_pearson = 0.9673`, `2,335` BH-FDR cis links).
- **Deterministic candidate eRegulons (`I6`, ADR-069):** offline MEME v4 scanning
  (`500 bp` peak centers), GC-quintile CMH cistrome enrichment, depth partial
  correlations, finite-population `(expr, prevalence)` null with pre-collapse BH-FDR,
  weighted GSEA leading-edge pruning, and paralog motif-sharing flags; `17.99x`
  hypergeometric edge fold-enrichment (`p = 2.93e-15`) and `19.78x` in shared signed
  modules vs `scenicplus==1.0.0` on `hc11` in `~166 s`.
- **Non-circular Cross-Modal Evidence Graph & `IntegrationNarrator` (`I7`, ADR-070):**
  `cross_modal_triangulation.py` separates `convergent_paired_multiome`
  (`shared_upstream=["peak_to_gene"]`, `independent_contrast_support=False`) from
  `convergent_cross_contrast`, `associative_cis_candidate`, `discordant_signal`, and
  `caveated_ambiguous`, rendered with 4 dual PNG (`300 DPI`) + editable SVG figures.

## Headline 2: Opt-In `HypothesisAgent` & Cross-Modal Adapter (`S0`–`S10`) — ADR-057..062, ADR-072

ARIA adds an opt-in (`--hypotheses` / `enable_hypotheses`) mechanistic hypothesis
generation layer quarantined strictly behind `SPECULATIVE` (`promotable=False`,
`assert_no_speculative_promotion`):

- **4 deterministic gates:** entity/direction/ledger grounding (`grounding.py`),
  falsifiable experimental design (`gates.py`), non-causal hedging (`gates.py`), and
  structured confounder/caveat ownership (`devils_advocate.py`), gated upstream by a
  closed `VerificationReceipt` (`W-CLAIM` + `W-LEDGER`).
- **5 pure evidence adapters (`bulk_rna`, `scrna`, `scatac`, `bulk_atac`,
  `integration`):** `adapters/integration.py` (`S10`, `ADR-072`) emits bounded typed
  signals from `peak_to_gene`, `eregulons`, and active non-excluded `shared` `mofa`
  RNA loadings (`direction="na"` honoring latent sign indeterminacy), collapsing
  shared cis lineage via code-owned `independence_key` so co-occurring `scATAC`
  `peak2gene`, `integration` `peak_to_gene`, and `integration` `eregulons` never
  double-count as independent lines of evidence.

## Headline 3: Bulk ATAC Beta Lane & Preprint Governance Hardening

- **Bulk ATAC-seq (`beta` + `requires_ack`, ADR-054):** FASTQ alignment (`bwa-mem2`),
  measured QC, MACS3 peak calling, replicate consensus, `bedtools coverage` count
  matrices, replicate-gated DESeq2 DA, local JASPAR2024 motif over-representation
  (validated on ENCODE K562 vs GM12878), and descriptive TOBIAS footprinting.
- **C1–C4 public-claim boundary & capsule portability:** untrusted prompt fencing,
  typed `SemanticFact` negation/polarity verification, provenance-linked claim-to-ledger
  compiler, causal estimand licensing, portable reproducibility capsules (`RO-Crate` +
  SQLite snapshot + lockfile hashes), atomic GEO/SRA staging (`E6`), and Track B3
  frozen open-weight reagent ablation (`Qwen3.6-35B-A3B-UD-Q4_K_M`, `25/25` `3/3`
  SHA-256 repeatable configs).

## Scope guard & preprint freeze status

DE/DA estimation core untouched (`_run_deseq2`, pre-registered FDR,
`lfcThreshold`-in-Wald, apeGLM). `scATAC`, `bulk_ATAC`, and paired `integration`
remain `beta` + `requires_ack` (`production_promoted=False`). Unpaired/diagonal/bulk
cross-modal integration, ChIP-seq, CUT&RUN, CUT&TAG, and Hi-C remain blocked or
scaffolded. Existing tags are never moved; the immutable `preprint-v1` artifact tag
remains gated (`freeze_gate.tag_action_authorized=False`) until the three independent
human-gold adjudication lanes complete.
