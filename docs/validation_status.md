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
| Chromatin scATAC alpha lane | Done for alpha | The scATAC matrix path is dispatchable behind explicit acknowledgement: measured QC, TF-IDF/LSI clustering, per-cluster DA, replicate-gated pseudobulk DA, local motif enrichment when resources exist, and chromatin narrative blocks. The modality remains alpha, not publication-grade autonomy |
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

## Alpha

Alpha modalities can dispatch only behind explicit acknowledgement and remain
review-required. They must surface missing resources, low replication, and
skipped lanes rather than fabricating output.

| Area | Status | Required before stable |
|---|---|---|
| scATAC-seq matrix workflow | Alpha + requires acknowledgement | Expert review of biological conclusions; broader fixtures and datasets; P2 regulatory layers are now input-gated alpha outputs (motif activity, gene scores, peak-to-gene, label-transfer hypotheses; Tn5 footprinting requires fragments + bias model); stable promotion requires independent review beyond the current HC11/synthetic/multi-sample validation evidence |

## Scaffolded / Roadmap

Scaffolded modalities are **dispatch-gated**: they cannot silently produce
publication-looking output. Direct calls to planned-but-absent scripts return a
structured `script_not_implemented`, and the integration scaffold returns an
explicit `NotImplemented` rather than any fabricated result.

| Area | Status | Required before stable |
|---|---|---|
| Bulk ATAC-seq | Scaffolded | Peak count matrix, DA, QC summaries, report section |
| ChIP-seq / CUT&RUN / CUT&TAG | Scaffolded | Clear assay-specific QC and peak interpretation |
| Hi-C / Micro-C | Scaffolded — dispatch OFF | Runs only under `ARIA_ALLOW_EXPERIMENTAL_HIC=1` and are stamped not-publication-grade; needs E2E validation + memory-safe fixtures |
| WNN / MOFA+ / peak-to-gene | Scaffolded | Stable standalone RNA + ATAC paths first; currently an explicit NotImplemented blocker, no fabricated weights or clusters |

Chromatin QC (`chromatin_qc.py`) emits only measured metrics. TSS enrichment and
FRiP are real when their inputs/resources exist and otherwise remain null with a
concrete skip reason, never a fabricated placeholder.

The v4.6 scATAC lane has been exercised on the local HC11 validation input
(`hc11_paired.h5mu`): ARIA reads ATAC modality `atac`, reports real dimensions
of 3,143 cells x 60,990 peaks, runs TF-IDF/LSI clustering, performs honest
single-sample DA where possible, and returns pseudobulk DA as skipped when
replicates/condition metadata are absent. Synthetic and real multi-sample
validation cover the pseudobulk DA execution path, but scATAC remains alpha.

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
