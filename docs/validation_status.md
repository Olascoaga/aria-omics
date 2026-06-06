# Validation Status

Last updated: June 6, 2026.

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

Current RNA/reporting baseline: `v4.5.4`, plus an ongoing post-`v4.5.4`
reliability / governance / reproducibility hardening pass (see below). scRNA
pseudobulk DE defaults to per-cluster FDR for the primary significance call while
still reporting global FDR for audit; significant-gene counts can differ from
pre-`v4.5.4` reports by design.

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
| Trajectory: PAGA + DPT | Beta | Validated on hippocampus subset; exploratory, not causal |
| Cell-cell communication: LIANA | Beta | Validated on GSE278576 annotated h5ad; `n_perms=1000` default for stable ranks |
| GEO/SRA connector | Beta | GSE183948 path validated; multi-organism (spike-in) organism inference added; public metadata remains heterogeneous |

## Scaffolded / Roadmap

Scaffolded modalities are **dispatch-gated**: they cannot silently produce
publication-looking output. Direct calls to planned-but-absent scripts return a
structured `script_not_implemented`, and the integration scaffold returns an
explicit `NotImplemented` rather than any fabricated result.

| Area | Status | Required before stable |
|---|---|---|
| scATAC-seq | Scaffolded (next, v4.6) | `.h5mu` entry path validated on the planned HC11 paired RNA+ATAC input; still needs LSI, clustering, differential accessibility, motifs, report section, fixtures |
| Bulk ATAC-seq | Scaffolded | Peak count matrix, DA, QC summaries, report section |
| ChIP-seq / CUT&RUN / CUT&TAG | Scaffolded | Clear assay-specific QC and peak interpretation |
| Hi-C / Micro-C | Scaffolded — dispatch OFF | Runs only under `ARIA_ALLOW_EXPERIMENTAL_HIC=1` and are stamped not-publication-grade; needs E2E validation + memory-safe fixtures |
| WNN / MOFA+ / peak-to-gene | Scaffolded | Stable standalone RNA + ATAC paths first; currently an explicit NotImplemented blocker, no fabricated weights or clusters |

Chromatin QC (`chromatin_qc.py`) already emits only measured metrics: TSS
enrichment and FRiP return null until a reference annotation / called peaks exist,
never a fabricated placeholder.

The v4.6 `.h5mu` entry path has been exercised on the planned local validation
input (`hc11_paired.h5mu`): ARIA reads ATAC modality `atac` and reports real
matrix dimensions of 3,143 cells x 60,990 peaks. This is not a production
scATAC validation. FRiP, TSS enrichment, complete QC, and pass/fail remain
explicitly not computed until the chromatin stack, peak calling, and
reference-backed TSS computation are implemented.

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
- `rna_de_per_cluster.py` is treated as optional on atlas-scale inputs. If it times
  out, the pipeline may still be scientifically valid when donor-level pseudobulk,
  pathway, communication, and trajectory outputs complete.

Latest reviewed hippocampus rerun produced real output tables and figures; newer
reports state that Leiden was skipped when input annotations are reused, and carry
the full provenance and evidence-tier manifest described above.
