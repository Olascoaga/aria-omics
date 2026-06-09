# ARIA v4.5 Benchmarking Protocol Freeze

---
status: active
source_of_truth_for: v45_benchmarking_protocol
last_updated: 2026-06-08
decision: ADR-030
---

## Position

The central contribution of ARIA is not a new differential expression test, but
a reproducible evidence-governance layer for agentic omics: one that preserves
established statistical behavior while reducing unsupported biological
interpretation.

This protocol freezes the v4.5 benchmarking spine. It is the methods backbone
for an RNA/preprint validation lane and closes the v4.5 product phase as a
publishable RNA + governance baseline. It does not block v4.6 scATAC
implementation unless the explicit goal changes to submitting the RNA preprint
before ATAC work resumes.

## Rationale

Recent agent benchmarks increasingly evaluate process, not only final answers,
because correct outcomes can arise from memorization, reward hacking, or flawed
reasoning. ARIA's evidence ledger, claim verification, refusal behavior, and
causal-language controls are aligned with this direction.

External benchmark use is supplemental. ARIA should run only a small compatible
subset of community tasks, primarily to show that its governance layer can
operate under independent task definitions. It should not attempt to win or
cover a full external benchmark.

Reference framing:

- BixBench: use as a supplemental external benchmark subset. Use the paper-safe
  description in main text: over 50 real-world computational-biology scenarios
  and nearly 300 associated open-answer questions. If using a released 205-question
  MCQ subset, cite the exact dataset/version in supplement.
- BiomniBench / process-level biomedical-agent benchmarks: cite for the
  motivation that outcome-only scoring misses flawed trajectories, method
  selection failures, and unsupported interpretation.

## Scope

The v4.5 benchmark answers three questions:

1. Does ARIA preserve standard RNA statistical behavior?
2. Does ARIA avoid statistically indefensible execution?
3. Does ARIA reduce unsupported or causally inflated biological narrative?

The benchmark is intentionally minimal. It must not become a second ARIA inside
ARIA.

## Main Figures

| Figure | Track | Main message |
|---|---|---|
| Fig 1 | A1 Bulk DE | ARIA preserves calibration and concordance against established bulk RNA references. |
| Fig 2 | A2 Pseudobulk scRNA | ARIA uses donor-aware pseudobulk analysis and does not treat cells as independent biological replicates. |
| Fig 3 | B1 DesignAgent | ARIA knows when to infer and when to block or escalate. |
| Fig 4 | B2 Claim/Narrative | ARIA reduces false narrative and causal overreach versus an ungoverned LLM baseline. |
| Fig 5 | B4 Null narrative | ARIA does not invent biology on null/permuted controls. |

Supplemental tracks: A3 composition, A4 annotation, A5 ORA/GSEA, A6
clustering, B3 multi-LLM invariance, and a small external BixBench/scBench
subset.

## Benchmark A: Statistical Validity

### A1 Bulk DE

Compare ARIA's bulk RNA path against established reference methods without
claiming identity or superiority. Report four axes separately:

1. FDR calibration on synthetic truth and ERCC/null controls.
2. LFC concordance, using Spearman/Pearson.
3. Ranking concordance, using top-k overlap, RBO, and Jaccard.
4. Significant-call concordance at FDR 0.05.

Expected result: high concordance, not identity. pyDESeq2 is a reimplementation
of DESeq2, and edgeR-QLF/limma-voom are related but not identical statistical
models.

SEQC/MAQC qPCR/TaqMan and ERCC spike-ins are external validation references,
not absolute truth.

SEQC/MAQC reference lane (executed): `aria/benchmarks/reference_seqc.py` +
`scripts/run_a1_seqc_maqc_benchmark.py` validate ARIA's real bulk DE path
against external TaqMan qPCR truth. The bundle is bootstrapped once from the
`seqc` Bioconductor data package by `scripts/fetch_seqc_maqc_reference.py`
(counts = a SEQC site's RefSeq gene table; truth = TaqMan log2(A/B) over ~1000
genes; kept out of the repo under `~/.aria/benchmarks/`). On the BGI RefSeq
counts (samples A/B, 5 reps each) ARIA scored: LFC concordance vs TaqMan
**Pearson 0.944 / Spearman 0.938** (830 genes), TaqMan-DE detection **AUC
0.893**, and titration monotonicity **97.5%** of TaqMan-DE genes ordered across
A→C→D→B. The lane is data-gated (`ARIA_SEQC_MAQC_BUNDLE`) and skips honestly when
the bundle is absent; nothing is fabricated (ADR-036).

Cross-site reproducibility (executed, ADR-037): `run_seqc_maqc_multisite`
(`scripts/run_a1_seqc_multisite_benchmark.py`) runs ARIA's A-vs-B DE at each of
the five Illumina SEQC sites (BGI, CNL, MAY, AGR, NVS) and reports the pairwise
log2FC concordance between sites — the SEQC reproducibility metric — plus each
site's TaqMan concordance. Executed result: mean off-diagonal cross-site Pearson
**0.980** (min 0.978) over all 10 site pairs (median 23,022 genes/pair), and a
near-constant per-site TaqMan Pearson **0.940–0.944**. ARIA's DE result is
effectively independent of the sequencing site. Artifact:
`docs/benchmark_results/a1_seqc_maqc/a1_seqc_multisite_v4.5.5.json`.

v4.5.5 executable artifact: the preliminary synthetic-truth A1 lane is
implemented by `scripts/run_a1_bulk_de_benchmark.py` and writes a versioned
manifest plus Fig 1 SVG under `docs/benchmark_results/`. It runs ARIA's real
bulk DESeq2 path with apeGLM enabled and reports the same four frozen axes.
External DESeq2/edgeR/limma comparator execution remains assigned to
`aria-bench-env`; do not cite the preliminary artifact as a superiority or
identity claim.

A1 also reports a permanent `lfc_threshold_frontier` axis (descriptive, not a
pass/fail gate): the same dataset is re-run through ARIA's real bulk DE path at
several Wald `lfcThreshold` values. `lfc_threshold=0` is the matched-null
DESeq2-equivalence reference (H0: LFC = 0) and is computed, not hardcoded;
higher thresholds test H0: `|LFC| <= thr`, trading recall for precision. On the
seed-11 synthetic truth ARIA recovers recall 0.808 / empirical FDR 0.085 at
`lfc_threshold=0` — matching the external DESeq2 comparator exactly — and
recall 0.525 / FDR 0.000 at the default `lfcThreshold=0.5` effect-size policy.
This isolates the recall difference versus DESeq2/edgeR/limma as ARIA's
deliberate, user-controlled policy rather than an engine difference; cite the
frontier, not a single conservative recall number (ADR-035).

Post-v4.5.5 comparator execution: `scripts/run_a1_external_comparators.py`
dispatches `aria/scripts/benchmark_a1_external_comparators.py` through
`EnvironmentManager` `stack="benchmark"` / `aria-bench-env`. The IPC runner
exports the same neutral A1 synthetic matrix and calls an R comparator script for
DESeq2, edgeR-QLF, and limma-voom, then scores each table against the known
truth. Local live execution wrote
`docs/benchmark_results/a1_external/a1_external_comparators_v4.5.5.json` with
all three methods complete. This is the synthetic-truth external-comparator
execution path; it is not a SEQC/MAQC/ERCC reference-data completion.

### A2 Pseudobulk scRNA

Use Kang et al. 2018 PBMC lupus control versus IFN-beta as the main scRNA
statistical figure. Compare ARIA donor-aware pseudobulk against muscat as the
multi-sample reference and include cell-level MAST/Wilcoxon-style baselines as
anti-pattern comparators.

Figure message: cell-level DE methods may inflate evidence when cells are
treated as independent replicates; ARIA preserves the donor/sample as the
inferential unit.

v4.5.5 executable artifact: the preliminary donor-aware A2 lane is implemented
by `scripts/run_a2_pseudobulk_benchmark.py` and writes a versioned manifest plus
Fig 2 SVG under `docs/benchmark_results/`. It validates ARIA's real pseudobulk
DE path on synthetic truth and runs a donor-heterogeneity null where a naive
cell-level Welch test treats cells as independent replicates. Kang + muscat
remains the external reference lane and requires local benchmark data plus
`aria-bench-env`.

### A5 ORA/GSEA

Keep three analyses distinct:

- ORA: hypergeometric/Fisher-style enrichment with local versus genomic
  background.
- Ranked GSEA: fgsea or decoupleR-style ranking.
- Sample-level scoring: GSVA only when the design supports it.

Do not compare these as if they answer the same estimand.

### A3/A4/A6 Supplemental

Use these only as sanity and rigor checks:

- A3 composition: propeller/scCODA/sccomp versus naive abundance tests.
- A4 annotation: macro-F1 versus reference labels when labels are defensible.
- A6 clustering: ARI/NMI and bootstrap stability when label/taxonomy limits are
  disclosed.

They are not preprint gates.

## Benchmark B: Governance

### B1 DesignAgent

Report three components separately, with a composite only as a summary:

- Correct inference rate.
- Correct refusal/block/escalation rate.
- Unsafe execution rate.

Unsafe execution rate is a headline metric and should be approximately zero.
It is defined as running an inferential analysis when the design is
statistically indefensible.

Adversarial cases include:

- fewer than three biological replicates per group;
- batch perfectly confounded with condition;
- ambiguous reference/control group;
- labels such as A/B with no biological semantics;
- denominator chosen alphabetically, for example Aged versus Young;
- continuous covariate disguised as a categorical group.

Report a confusion matrix by decision type.

### B2 Claim/Narrative

Operational unit:

> A scientific claim is any sentence that asserts a biological, statistical,
> methodological, or interpretive conclusion about the analyzed dataset.

Excluded from the denominator:

- generic background/context;
- method descriptions without interpretation;
- figure-navigation text;
- administrative report text.

Claim segmentation is frozen and double-annotated. Agreement is measured for
segmentation and labeling, not only final labels.

Claim schema:

```text
claim_id
sentence
subject
predicate
object
evidence_ids
evidence_tier
claim_type: descriptive | comparative | mechanistic | causal | speculative
label_expert_1
label_expert_2
adjudicated_label
```

Failure taxonomy:

- unsupported;
- overclaim;
- fabricated;
- causal inflation;
- missing caveat.

Primary metrics:

```text
False narrative rate =
  (unsupported + fabricated + overclaim) / total scientific claims

Causal overreach rate =
  causal_or_mechanistic_claims_without_matching_evidence / total causal_or_mechanistic_claims
```

Main ablation arms:

1. ARIA governed.
2. ARIA without Claim Compiler / evidence guards.
3. Naive LLM report.
4. Template-only report.

Expected message: ARIA governed approaches the safety of deterministic
templates while preserving richer synthesis and caveat integration, and remains
well below the naive-LLM false-narrative rate.

The gold standard is human annotation with adjudication. Any automated scorer
must be independent of ARIA's Claim Compiler and reported against human labels.
At least one dataset should be recent/private enough to reduce training-data
contamination risk.

### B3 Multi-LLM Invariance

Supplemental unless it is clean and cheap:

- 3 datasets;
- 3 models;
- 3 repetitions;
- temperature 0;
- frozen input manifests.

Compare hashes of statistical outputs, DE tables, claim counts/support, and
non-scientific wording variation.

### B4 Null Narrative

Use label permutations and null/spike-in controls tied to W-CALIB. Measure
whether ARIA reports no supported signal rather than inventing a biological
story.

Primary metric: fabricated or unsupported narrative rate on null controls.

## Frozen Implementation Order

1. Benchmark manifest schema, compatible with the W-CALIB badge manifest.
2. A1 synthetic bulk DE and preliminary Fig 1. **Done in v4.5.5 for the
   ARIA-path synthetic-truth lane.**
3. `aria-bench-env.yml` for external comparators. Keep R/benchmark packages out
   of production RNA/chromatin environments and call them through JSON IPC.
   **Env + A1 IPC runner scaffold implemented; live R comparator execution
   requires local `aria-bench-env`.**
4. A2 Kang + muscat. **Preliminary ARIA-path donor-aware lane done in v4.5.5;
   external Kang + muscat remains pending.**
5. B1 adversarial design corpus, about 30 cases, no heavy downloads.
6. B2 claim schema and manual scorer CSV/rubric.
7. B4 null narratives.
8. Figures and tables from the start, not as an afterthought.
9. Small external BixBench/scBench subset, 2-3 compatible RNA/scRNA tasks.
10. B3 multi-LLM invariance only if it does not delay the core paper.

B2 and B4 outrank B3.

## Release Gate

The v4.5 product line is closed by freezing this benchmark protocol and the
RNA/governance validation story. Full execution of every benchmark track is a
preprint milestone, not a prerequisite for starting v4.6 scATAC.

If the immediate product goal changes to "submit the RNA preprint before
scATAC", then A1, A2, B1, B2, and B4 become the minimum execution gate.

## Guardrails

- Agreement is not truth; disclose proxy endpoints.
- Do not over-sell Benchmark A. Concordance with DESeq2-like behavior validates
  non-regression, not superiority.
- Penalize correct numbers for the wrong reasons via evidence/process metrics.
- Keep the external benchmark small and supplemental.
- Do not condition the preprint on v4.7 integration.

## Product Phases

| Phase | Tracks | Preprint gate |
|---|---|---|
| RNA / v4.5 | A1, A2, A3, A5, B1, B2, B3, B4, external subset | Yes, if submitting RNA preprint now |
| v4.6 ATAC | scATAC DA/QC/LSI/motifs validation | Revised manuscript / next validation lane |
| v4.7 integration | WNN/MOFA+/peak2gene and multi-omic B tracks | Follow-up |
