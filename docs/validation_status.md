# Validation Status

Last updated: June 12, 2026.

ARIA uses explicit validation boundaries so users can distinguish mature
workflows from beta analysis paths and implementation scaffolds.

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
| scRNA single-sample | Validated | PBMC 3k and GSE278576 single-sample runs |
| scRNA multi-sample | Validated | GSE278576 3-donor concat + Harmony workflow |
| Processed h5ad pseudobulk | Validated | 40-donor hippocampus rerun; report review confirmed pseudobulk, ORA, LIANA, trajectory, figures, and TSV exports |
| Integrated RNA biological discussion | Validated for RNA-only reports | Bulk RNA and scRNA synthesis blocks are generated from structured ARIA outputs, evidence-verified, claim-tiered, and rendered through the normal NarrativeAgent governance |

Current RNA/reporting baseline: release tag `v4.6`, with post-`v4.6` hardening
on `main` through the June 12, 2026 scRNA audits and production E2E closure. The
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
| scRNA-lane production audit | Done | Closed B-PB1, B-DD1, B-QC1/B-QC2, A-MARK1, A-CMT1/A-CLUST1, and B-TRAJ1. Production pseudobulk uses raw QC counts, marker claims remain descriptive, QC is data-intrinsic, and DPT is skipped with `root_unresolved` when no defensible root exists |
| scRNA annotation/integration audit | Done | Closed N-ANNO1/N-ANNO2/N-ANNO3, N-QC1, and N-INT1. CellTypist confidence is genuine, low-confidence cell-type labels cap pseudobulk block confidence, default immune model fallback is disclosed, count-MAD QC is log-space, and destructive integration overcorrection is blocking |
| scRNA production E2E verification | Done | Real raw 10X pbmc3k runs verified the above report surfaces. Raw 10X MEX directories are collapsed to one sample, QC failures render as errors, and multi-donor pseudobulk verified `count_source=raw_counts` with the handoff |

## Reliability, Governance & Reproducibility

These cross-cutting guarantees back every validated path and are enforced in CI:

| Guarantee | Status | How |
|---|---|---|
| Deterministic narrative | Done | LLM calls at `temperature=0` + fixed seed; model / tier / cache hit / degraded fallback recorded per report |
| Adversarial review of claims | Done | Deterministic, LLM-free devil's advocate enumerates the standard confounders (batch, ambient RNA, doublets, composition shift, low replication) per claim |
| Planned-vs-run ledger | Done | Report reconciles planned vs executed analyses; a partial run is visible, not silent |
| Single-source provenance | Done | version + git commit + dirty state + input SHA-256 + per-stage parameter hashes + dependency lockfiles + LLM usage |
| Design honesty | Done | confirmed covariates in the fitted DESeq2 formula; no alphabetical reference contrast; no filename-fallback design in production; user thresholds propagated end-to-end |
| No fabricated science | Done | repo-wide anti-fabrication guard (no placeholder matrices, ungated mock "successes", or hash-derived metrics) |
| No hardcoded biology in runtime fallbacks | Done | biological marker panels, ligand-receptor fallback tables, tissue keyword maps, and mechanism labels from feature names/signs are not used as production evidence |
| Evidence-linked claims | Done | report claims pass strict evidence verification, causal-language guards, run-ledger linkage, and deterministic devil's-advocate review before rendering |
| Typed IPC contracts | Done | every dispatchable script validates inputs/outputs before and after the subprocess |
| Blocking CI (3 tiers) | Done | PR (guards + unit + contracts) / main (real pyDESeq2 recovery benchmark) / release (Docker env solve + in-image benchmark) |
| Air-gapped mode | Done | `ARIA_AIR_GAPPED=1` governs **all** egress (LLM + pathway ORA + GEO/SRA connectors), not just the LLM, and redacts failed-run input/error archives |
| Sensitivity checkpoint | Done | inputs are classified for clinical/PHI-like fields and quasi-identifiers before CP1; the user is always offered an air-gapped opt-in (never auto-disabled) |
| Local versioned enrichment | Done | over-representation analysis runs locally by default against versioned GMTs (library + release + SHA-256 recorded); Enrichr is opt-in and skipped honestly when egress is blocked |
| Deterministic packaging | Done | PEP 621 `pyproject.toml` with a single version source and per-platform lockfiles; `requirements.lock` core fallback |
| Hermetic build + supply chain | Done | per-modality Docker images with a `.dockerignore` (no `.git`/private `memory/`/caches in layers), gitleaks secret scan, and a CycloneDX SBOM of the RNA image; image digest stamped in the report |
| Secret hygiene | Done | installer reads API keys without echo; `aria doctor --secrets`/`--llm` classify/mask keys and flag committed credentials (no LLM call) |
| Per-modality readiness gate | Done | an assay capability matrix marks each modality green/yellow/red; green auto-dispatches, yellow needs explicit acknowledgement, red is removed from dispatch |

## Validated / Beta

| Area | Status | Notes |
|---|---|---|
| Bulk RNA FASTQ preprocessing | Beta | Scripted path exists; dependency and real-data coverage should expand |
| Trajectory: PAGA + root-gated DPT | Beta | PAGA validated on hippocampus subset; DPT requires precomputed `iroot`, an explicit matching root label, or a generic progenitor/stem/precursor label. If no root is available, pseudotime is skipped with `root_unresolved`; trajectory remains exploratory, not causal |
| Cell-cell communication: LIANA | Beta | Validated on GSE278576 annotated h5ad; `n_perms=1000` default for stable ranks |
| GEO/SRA connector | Beta | GSE183948 path validated; multi-organism (spike-in) organism inference added; public metadata remains heterogeneous |
| scATAC-seq matrix workflow (de-alpha v4.7.0 / ADR-048) | Beta + requires acknowledgement | Scoped beta: QC/clustering, motif, and the replicate-gated pseudobulk condition-DA lane are externally concordant (SnapATAC2 on HC11; edgeR-QLF/limma-voom/R-DESeq2 on 5v5 GSE278576 donors — ARIA ⊆ R-DESeq2 recall 1.000, LFC Spearman 0.61-0.76). The per-cluster Wilcoxon marker path stays caveated as single-sample-fragile. **Peak-to-gene link recovery is promoted to beta (P4.2 / ADR-050):** externally concordant with Signac LinkPeaks on HC11 (GSE278576, 3143 paired cells) — Spearman 0.62 and 99.98% sign agreement on the 14,046 shared links; moderate set overlap (Jaccard 0.20) is method/threshold-driven (raw-count Pearson vs TF-IDF + permutation background), links stay associative. The remaining P2 regulatory layers (motif activity, gene scores, label-transfer hypotheses) stay input-gated scaffold; gene activity needs fragments and is deferred to the fragment track (P4.3); Tn5 footprinting refuses uncorrected output (honest not-run). Stable promotion still needs independent expert review of biological conclusions |

## Alpha

No modalities are currently at alpha: scATAC was promoted to scoped **beta** in
v4.7.0 (ADR-048). Any future alpha modality dispatches only behind explicit
acknowledgement, remains review-required, and must surface missing resources, low
replication, and skipped lanes rather than fabricating output.

## Scaffolded / Roadmap

Scaffolded modalities are **dispatch-gated**: they cannot silently produce
publication-looking output. Direct calls to planned-but-absent scripts return a
structured `script_not_implemented`, and the integration scaffold returns an
explicit `NotImplemented` rather than any fabricated result.

| Area | Status | Required before stable |
|---|---|---|
| Bulk ATAC-seq | Beta slice | Measured QC + MACS3 peak calling are acknowledgement-gated beta; peak count matrix + DA remain scaffolded |
| ChIP-seq / CUT&RUN / CUT&TAG | Scaffolded | Clear assay-specific QC and peak interpretation |
| Hi-C / Micro-C | Scaffolded — dispatch OFF | Runs only under `ARIA_ALLOW_EXPERIMENTAL_HIC=1` and are stamped not-publication-grade; needs E2E validation + memory-safe fixtures |
| WNN / MOFA+ / peak-to-gene | Scaffolded | Stable standalone RNA + ATAC paths first; currently an explicit NotImplemented blocker, no fabricated weights or clusters |

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
Wilcoxon marker path stays caveated as single-sample-fragile; scATAC is beta +
`requires_ack`, not autonomous publication-grade.

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
- `rna_de_per_cluster.py` is treated as optional on atlas-scale inputs. If it times
  out, the pipeline may still be scientifically valid when donor-level pseudobulk,
  pathway, communication, and trajectory outputs complete. DPT pseudotime may be
  absent by design when no defensible root is available.

Latest reviewed hippocampus rerun produced real output tables and figures; newer
reports state that Leiden was skipped when input annotations are reused, and carry
the full provenance and evidence-tier manifest described above. The integrated
RNA discussion is part of the governed report path; cross-modal RNA+ATAC
synthesis remains deferred until validated chromatin outputs exist.
