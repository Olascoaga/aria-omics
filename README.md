# ARIA

### Agentic Research Intelligence for -omics Analysis

> *You ask the biological question. ARIA does the rest.*

![Version](https://img.shields.io/badge/version-4.7.0-blue)
![Status](https://img.shields.io/badge/status-active-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.11-blue)

---

## What is ARIA?

ARIA is an open-source agentic system for supervised omics analysis. Instead
of only exposing pipeline commands, ARIA asks for the biological question,
confirms the experimental design, runs modality-specific code, records
decisions, and writes a report grounded in real output files.

**Current validation boundary:**

- **Validated (controlled + small real datasets):** bulk RNA-seq count
  matrices; scRNA-seq single-sample and multi-sample workflows; processed
  `.h5ad` pseudobulk workflows with design metadata in `obs`. "Validated" means
  these paths were exercised on controlled synthetic data and small real
  datasets with reviewed reports and a numerical ground-truth benchmark — it is
  **not** a claim of publication-grade results for a specific study; a domain
  expert must still review the design, the fitted model, and the conclusions
  before publication.
- **Validated / beta:** bulk RNA-seq FASTQ preprocessing, trajectory
  summaries (PAGA plus DPT only when a defensible root is available), LIANA
  cell-cell communication, processed `.h5ad` recovery, and GEO/SRA connector
  paths.
- **Beta (requires acknowledgement):** scATAC-seq matrix workflow — QC/clustering,
  motif enrichment, and the replicate-gated pseudobulk condition-DA lane are
  beta-grade (de-alpha in v4.7.0 / ADR-048, externally concordant with SnapATAC2
  and edgeR/limma/DESeq2); the per-cluster Wilcoxon marker path stays caveated as
  single-sample-fragile. Bulk ATAC-seq is open as a V47 beta slice for measured
  QC + MACS3 peak calling, and comparison requests can now build a scaffolded
  peak-by-sample count matrix; bulk ATAC differential accessibility remains
  scaffolded. Both lanes are dispatch-gated behind explicit acknowledgement.
- **Scaffolded / roadmap:** ChIP-seq, CUT&RUN / CUT&TAG, full Hi-C / Micro-C
  workflows, and multimodal WNN/MOFA+ integration.

**ARIA produces:**

- Pre-analysis quality audit with actionable warnings before expensive compute
- Differential expression with explicitly confirmed contrasts (never an
  alphabetically guessed reference), covariate/batch-adjusted DESeq2 designs
- Pathway enrichment (ORA + GSEA) per contrast, with explicit gene-set background
- Publication-ready figures: volcano plots, PCA/MDS, ORA dotplots, GSEA running sums
- HTML report in paper/publication style with embedded figures
- Manuscript-ready methods section with the exact fitted design and parameters
- Evidence-tiered claims with a per-claim manifest and an adversarial caveat pass
- A planned-vs-run ledger so a partial run is visible, not silent
- Full provenance: version + git commit + dirty state + input SHA-256 + per-stage
  parameter hashes + dependency lockfiles + LLM usage, in HTML and
  `methodology.json`
- Reproducible decision log (every threshold choice stored in SQLite)

---

## Documentation

Start with [docs/README.md](docs/README.md).

Core public docs:

- [Architecture overview](docs/architecture/overview.md)
- [Design principles](docs/architecture/design_principles.md)
- [Reporting and outputs](docs/architecture/reporting_and_outputs.md)
- [Validation status](docs/validation_status.md)
- [Installation guide](docs/INSTALLATION.md)

Workflows:

- [Bulk RNA-seq](docs/workflows/bulk_rna.md)
- [Single-cell RNA-seq](docs/workflows/scrna.md)
- [Pseudobulk scRNA-seq from h5ad obs metadata](docs/workflows/pseudobulk_scrna.md)
- [Trajectory analysis](docs/workflows/trajectory.md)
- [Cell-cell communication](docs/workflows/cell_communication.md)

Roadmaps:

- [Chromatin workflows](docs/workflows/chromatin_roadmap.md)
- [Hi-C / Micro-C workflows](docs/workflows/hic_roadmap.md)
- [Multimodal integration](docs/workflows/integration_roadmap.md)

Diagrams are stored as Mermaid files under [docs/diagrams](docs/diagrams/).

---

## Current Status — June 2026

```
INFRASTRUCTURE                     STATUS
────────────────────────────────────────────────────
MessageBus (durable, per-run replay) done
ARIAMemory (SQLite hierarchical)   done
LLMProvider (deterministic+airgap) done — temp=0, fixed seed, local-only mode
ParameterAdvisor (3-layer)         done
EnvironmentManager (typed IPC/JSON) done — contract for every dispatchable script
Anti-fabrication guard (CI gate)   done — no placeholder/mock/hash-derived output
Provenance stamping (single source) done — version + commit + dirty + workflow hash
DebateCouncil + devil's advocate   done — deterministic adversarial review

AGENTS                             STATUS
────────────────────────────────────────────────────
OrchestratorAgent                  done
DataAuditAgent                     done
DesignAgent (interactive, v4.0)    done
AuditAgent (quality linter, v4.1)  done
SetupAgent                         done
BulkRNAAgent (DESeq2)              done  ✓ H9 BMAL1/REV-ERBα
  rna_bulk_de.py                   done
  rna_pathway_viz.py               done
scRNAAgent                         done  ✓ PBMC 3k + GSE278576 multi-donor
  rna_qc.py  (MAD + Scrublet)      done — per-sample mode with sample_id
  rna_concat.py (multi-sample)     done — inner-join concat, raw preserved
  rna_integration.py (Harmony)     done  ✓ validated on 3 donors
  rna_advise_resolution.py         done
  rna_clustering.py                done — Leiden or reused obs annotations
  rna_celltypist.py                done
  rna_de_per_cluster.py            optional — may time out on large atlases
  rna_pseudobulk_de.py             done — between-condition DE via pyDESeq2
  rna_pathway_per_cluster.py       done — also used for per (group, comp) ORA
  rna_trajectory.py (PAGA + root-gated DPT)
                                      beta — validated on hippocampus subset
  rna_cellcomm.py (LIANA)          beta — validated on GSE278576
ChromatinAgent                     beta — scATAC and bulk ATAC dispatch require explicit ack
  chromatin_qc.py                  done — measured-only QC (no fabricated TSS/FRiP)
  chromatin_lsi_clustering.py      beta — TF-IDF/LSI clustering (SnapATAC2-concordant)
  chromatin_diffacc.py             beta — replicate-gated pseudobulk DA (edgeR/limma/DESeq2-concordant); per-cluster Wilcoxon markers caveated
  chromatin_motifs.py              beta — local motif enrichment when resources exist
  chromatin_peaks.py (MACS3)       beta for bulk ATAC peak calling; other assay use remains scaffolded
  chromatin_peak_counts.py         scaffold — bulk ATAC peak x sample matrix; DA still gated
GenomeArchAgent                    scaffolded — dispatch OFF by default
  hic_inspect.py                   done — needs ARIA_ALLOW_EXPERIMENTAL_HIC=1
  hic_qc_and_balance.py            done — runs are stamped not-publication-grade
  hic_topology.py (out-of-core)    scaffolded
NarrativeAgent (HTML report)       done  ✓ paper theme
  Narrative kernel (evidence cards) done — validated NarrativeBlock objects
  Claim Compiler (evidence tiers)   done — per-claim tier + manifest (X14)
  Run ledger (planned vs executed)  done — partial runs surfaced, not silent
  Devil's advocate (deterministic)  done — confounder check per claim
IntegrationAgent (WNN + MOFA+)     scaffolded — dispatch-gated, emits no
                                   fabricated output (explicit NotImplemented)
GEO/SRA connectors                 done   ✓ GSE183948 validated
```

The current release baseline is `v4.7.0`: bulk RNA + scRNA core paths are closed
for practical use, publication-readiness provenance is embedded in reports,
raw-ingestion planning/conversion is available for supported 10X inputs,
reports are composed from validated modality blocks for scRNA and bulk RNA, the
interactive "ARIA Control Center" (Textual TUI) is opt-in over the canonical
headless path, and the scATAC matrix workflow is a scoped **beta** path behind
explicit acknowledgement (de-alpha in v4.7.0 / ADR-048).

`v4.5.3` tagged the pre-ATAC integrity freeze: centralized version metadata,
installer secret hygiene, registry-integrity checks, scaffold dispatch gating,
typed script IPC contracts, design-matrix validation before DESeq2, synthetic
ground-truth DE benchmarking, scientific QC red flags, and the **Claim
Compiler**. `v4.5.4` adds scientific-honesty hardening: pseudobulk scRNA now
defaults to per-cluster FDR for the primary significance call while still
reporting global FDR as a secondary audit diagnostic. Power is reported against
the decision rule actually in force — under the per-cluster default, against
each block's effective per-cluster-family alpha (the whole-experiment global-BH
alpha is kept only as a secondary diagnostic) — and log-normalized count
recovery is visibly low-confidence. Pseudobulk significant-gene counts can
differ from pre-`v4.5.4` reports by design.

`v4.5.5` adds the first executable artifact from the frozen v4.5 benchmarking
plan: preliminary Benchmark A1 bulk-DE validation against synthetic truth, with
FDR calibration, LFC concordance, ranking concordance, and significant-call
concordance reported in a versioned manifest and simple Fig. 1 SVG. External R
comparators remain in the separate `aria-bench-env` lane.

Since `v4.5.4`, ARIA has gone through a focused reliability, governance, and
reproducibility hardening pass on the validated RNA baseline, followed by the
`v4.6` scATAC alpha line and post-`v4.6` scRNA production fixes. Shipped so far:

- **Deterministic, auditable narrative** — every LLM call runs at
  `temperature=0` with a fixed seed, and each report records which model/tier
  answered and whether it was a cache hit or a degraded fallback. A
  deterministic *devil's advocate* challenges every associative-or-stronger
  claim with the standard technical confounders (batch, ambient RNA, doublets,
  composition shift, low replication) and records which were addressed.
- **Planned-vs-run ledger** — every report reconciles the analyses the plan
  called for against the ones that actually executed, so a partial run is
  visible instead of silent.
- **Single-source provenance** — each report stamps the exact version, git
  commit, dirty state, and workflow hash, all derived from one version source.
- **Design honesty** — bulk DE now honors confirmed batch/covariates in the
  fitted DESeq2 formula, refuses to pick a reference contrast by alphabetical
  order, refuses to silently infer the design from file names, and propagates
  the user's confirmed significance thresholds end-to-end.
- **No fabricated science** — a repo-wide guard fails the build if any script
  returns placeholder matrices, ungated mock "successes", or hash-derived
  metrics; every dispatchable analysis script carries a typed IPC contract; and
  Hi-C / integration scaffolds cannot emit publication-looking output by
  default (Hi-C requires an explicit `ARIA_ALLOW_EXPERIMENTAL_HIC=1` opt-in and
  is stamped not-publication-grade).
- **Governed execution** — a blocking three-tier CI (PR / main / release) runs
  the guards, the real pyDESeq2 numerical-recovery benchmark, and a Docker
  env-solve; an air-gapped mode (`ARIA_AIR_GAPPED=1`) keeps sensitive runs
  local-only and redacts failed-run archives.
- **scRNA production-lane hardening** — production pseudobulk keeps raw
  QC-filtered counts through an explicit `counts_data_path` handoff, per-cluster
  marker discovery is reported as descriptive because it is cluster-defined and
  double-dipped, QC thresholds are data-intrinsic and use log-space MAD bounds
  for count metrics, CellTypist fallback models and low-confidence labels are
  disclosed in reports, overcorrection signatures can block integration quality,
  raw 10X MEX directories are treated as one sample, and failed QC blocks render
  as failures rather than success.

**End-to-end validated on bulk RNA-seq** — human H9 cells (3 conditions × 3 replicates):

- All pairwise contrasts generated: BMAL1 vs WT, REV-ERBα vs WT, BMAL1 vs REV-ERBα
- DESeq2 + ORA (GO_BP, KEGG, Reactome) + GSEA per contrast
- Publication-ready HTML report with embedded volcano plots and pathway dotplots

**End-to-end validated on scRNA-seq** (v4.3.3) — `tests/test_scrna_e2e.py`:

- **PBMC 3k** (single-sample) → `Immune_All_Low` CellTypist model: T helper,
  classical monocytes, CD16⁺ NK, B cells (4 clusters at Leiden 0.2). Per-cluster
  ORA recovers expected biology (B cells → MHC II antigen presentation;
  monocytes → neutrophil degranulation; NK → cytotoxicity).
- **GSE278576 hippocampus (hc11)** (single-sample) →
  `Adult_Human_PrefrontalCortex` model: OPC, Oligo, Astro, Micro, InN VIP
  (5 clusters at Leiden 0.2). Per-cluster ORA: glutamatergic synapse,
  axon development.
- **GSE278576 hippocampus (3 donors)** (multi-sample) — exercises
  `rna_concat` + `rna_integration` (Harmony): 11,783 cells × 22,406 shared
  genes; Harmony silhouette −0.047 → −0.067 (lower = better mixing); 9
  clusters at Leiden 0.2 (silhouette 0.677) annotated as Oligo, Endo, L2-3
  excitatory neurons, InN VIP, Astro AQP4, Astro GFAP, Micro, etc.
- **GSE278576 hippocampus (consolidated, 40 donors, 295k cells)** —
  pseudobulk DE on a Seurat-derived `.h5ad` (lognorm-counts recovered
  in-stream from `nCount_RNA`): age_group `80-100 vs 20-39` within each
  of 18 cell-type subclasses, Gender as covariate, pyDESeq2 design
  `~ age_group + Gender`. Recovered ORA signal matches the brain-aging
  literature: electron transport chain / oxidative phosphorylation
  dysfunction across Chandelier / LAMP5 / PVALB / SUB neurons and
  astrocytes; type I interferon signaling in VLMC; Parkinson disease
  pathway in SUB.
- **GSE278576 hippocampus full aging rerun** — validated on a 40-donor
  processed `.h5ad` with 295,033 starting cells and 242,405 cells retained
  after QC. ARIA reused the existing `obs['subclass']` annotations, skipped
  Leiden clustering, and treated donor-level pseudobulk DE as the primary
  inferential layer. The report exported UMAPs, pseudobulk DE tables, ORA
  tables, LIANA interactions, and trajectory summaries. PAGA can be present
  without DPT when no defensible pseudotime root is available. The legacy
  cell-level `rna_de_per_cluster.py` stage may time out on atlas-scale inputs;
  this is non-blocking when pseudobulk outputs are available.

**Publication-grade HTML report for scRNA / pseudobulk** (v4.3.5) — run
`tests/test_scrna_e2e.py --emit-html` and ARIA writes a self-contained
HTML to `<workspace>/report/` containing: UMAP coloured by cell-type and
condition, per-cell-type DE summary bar, per-(group×comparison) DE table
sorted by significance with top up/down genes, ORA dotplots for the top
cell types across GO_BP / KEGG / Reactome, Methods section with exact
thresholds and design formula, and the parameter-decisions log.

**Trajectory analysis (PAGA + root-gated DPT)** (v4.3.6) — preprocessed h5ads can
be dispatched with `--trajectory-h5ad PATH --trajectory-groupby COL
--trajectory-root TYPE`. PAGA connectivity can run without a root; DPT
pseudotime is computed only when ARIA has a defensible root (`iroot`, an
explicit matching `--trajectory-root`, or a generic progenitor/stem/precursor
label). When no root is available, ARIA reports `root_unresolved` instead of
choosing low-complexity cells as a proxy. Validated on a 50k-cell OPC + Oligo +
Astro subset from the hippocampus dataset: with OPC selected as root, DPT
pseudotime orders OPC (0.216) → Oligo (0.239) → Astro (0.367). PAGA reports
`max_connectivity = 0.00231` — characteristic of mature / non-developmental
tissue — and the HTML embeds both a normal and a log-scaled cluster graph so
weak adult-tissue edges remain visible without the report making false
developmental claims. RNA velocity is automatically skipped (the dataset lacks
spliced / unspliced layers); re-quantification via velocyto or kb-python `nac`
mode is documented in the Methods section as the path to enable it.

**Cell-cell communication (LIANA)** (v4.3.7) — dispatch with
`--cellcomm-h5ad PATH --cellcomm-groupby cell_type_celltypist`.
Validated on the GSE278576 multi-sample annotated h5ad (9 cell types):
LIANA `rank_aggregate` (now `n_perms=1000` by default for stable
publication-grade ranks). Autocrine pairs (source == target) are excluded a
priori (~1.5k of ~13.5k rows on this dataset). The script auto-falls back to `specificity_rank` when LIANA
emits all-NaN `magnitude_rank` (a recent-version quirk). Top non-
autocrine hits recovered are classical glia-neuron signaling axes:
APOE → TREM2 (Astro → Micro), C3 → NRP1 (Micro → L2-3 neurons),
VCAN → TLR2 (OPC → Micro), NRG1 → MS4A4A (InN VIP → Micro). HTML
includes a sender × receiver interaction-count heatmap, the top-N
L-R bar chart, and a sortable table with CellPhone p-values.

---

## Design Principles

**Local-first** — Your data never leaves your machine. An air-gapped mode
(`ARIA_AIR_GAPPED=1`) restricts the LLM layer to local models, refuses cloud
calls when no local model is configured, and redacts failed-run input archives.

**Language as interface** — Ask a biological question. ARIA translates it
into an analysis plan, executes it, and explains what it found.

**Supervised autonomy** — Checkpoints let you review and correct before
critical decisions. The DesignAgent walks you through experimental design
interactively (groups, organism, covariates, batch factors) before any
computation starts.

**Inspect before investing** — The AuditAgent runs a quality linter after
design confirmation and before expensive dispatch. It flags sample swaps,
outliers, batch dominance, and low alignment rates — with specific,
actionable recommendations.

**Honest uncertainty** — Every finding carries a confidence level
(HIGH / MEDIUM / LOW / INSUFFICIENT). ARIA tells you when data is
ambiguous or underpowered.

**Evidence-tiered claims** — A deterministic Claim Compiler classifies every
biological claim by the structured evidence that actually supports it —
*descriptive → associative → weak-mechanistic → strong-mechanistic →
causal-experimental* — and caps the language the report may use. Observational
omics is reported as association, not causation, unless the design is
interventional; claims whose wording exceeds their evidence tier are flagged.
Each claim ships with an evidence manifest in `methodology.json`.

**Reproducible by construction** — The narrative is deterministic
(`temperature=0`, fixed seed), every report stamps its exact version, git
commit, dirty state, input hashes, and dependency lockfiles, and a blocking CI
runs a real numerical-recovery benchmark plus an anti-fabrication guard on every
change. Nothing reaches a report that the code cannot reproduce.

**Institutional memory** — Every approved parameter decision is stored
in a local SQLite database. Over time, ARIA learns your lab's analytical
preferences and cites historical decisions in its justifications.

**Dependency isolation** — Each analytical stack runs in its own Conda
environment. IPC via JSON files prevents C-library conflicts between
scanpy, cooler, MACS3, and pysam.

**Token-efficient** — Inter-agent communication uses a compact internal wire
format. Only user-facing outputs are in normal prose.

**Provider-agnostic** — Works with Claude (Anthropic), Gemini (Google),
or local models via Ollama. Switch providers in one config line.

---

## Architecture

```
ARIA
  OrchestratorAgent       Parses question → plan → coordinates all agents
  DataAuditAgent          Auto-detects data types          [CHECKPOINT 1]
  DesignAgent             Interactive experimental design  [CHECKPOINTS 2.1–2.6]
    Groups → Organism → Factor → Batch → Pseudoreplication → Confirm
  Analysis Plan           User confirms or modifies plan   [CHECKPOINT 2]
  Threshold Tuning        Optional DE threshold profile     [CHECKPOINT 3]
  AuditAgent              Quality linter before dispatch   [CHECKPOINT 3.5 if blocking]
    · Replicate correlation (outlier / swap detection)
    · PCA batch dominance
    · STAR alignment rate
  SetupAgent              Environment/genome readiness check before dispatch
  BulkRNAAgent            DESeq2, all pairwise contrasts, ORA, GSEA
  scRNAAgent              QC, clustering, annotation, DE
  ChromatinAgent          scATAC matrix workflow [beta, requires ack]
                          bulk ATAC QC + MACS3 peaks [beta, requires ack]
                          bulk ATAC peak-count matrix [scaffold, DA gated]
                          ChIP + CUT&RUN + CUT&TAG [scaffolded]
  GenomeArchAgent         HiC, TADs, loops, compartments   [scaffolded]
  IntegrationAgent        Conditional multimodal synthesis  [scaffolded]
  NarrativeAgent          HTML report + methods section
  Final Review            Accept / revise / export          [CHECKPOINT 5]

  LLMProvider       Universal LLM abstraction (Anthropic / Gemini / Ollama)
  ContextManager    4-step degradation cascade for local models
  ParameterAdvisor  3-layer hyperparameter decisions + institutional memory
  EnvironmentManager IPC via JSON, isolated Conda stacks per modality
  Claim Compiler    Deterministic evidence-tiering + adversarial caveat pass
  DebateCouncil     Internal LLM peer review: Proposer vs Critic (2–3 rounds)
  ARIAMemory        Hierarchical SQLite: Wings / Halls / Rooms / Findings
  MessageBus        Durable per-run pub/sub with compact internal messages
  TUI               Terminal interface (Rich)
```

---

## Checkpoints

| # | When | What ARIA asks |
|---|------|----------------|
| 1 | After data scan | "This is what I found — is it correct?" |
| 2.1 | Design: groups | Confirm or correct detected sample groups |
| 2.2 | Design: organism | Select or type organism + genome assembly |
| 2.3 | Design: factor | Choose the main experimental variable |
| 2.4 | Design: batch | Declare batch covariates if present |
| 2.5 | Design: pseudoreplication | Confirm biological vs technical replicates |
| 2.6 | Design: confirm | Review full experimental design before planning |
| 2 | Analysis plan | "Here is my plan — shall I proceed?" |
| 3 | Parameter tuning | Select DE threshold profile (strict / standard / exploratory) |
| 3.5 | Quality audit | "Blocking issues found — proceed anyway or cancel?" |
| 5 | Final report | "Analysis complete — review or export?" |

---

## Installation

```bash
git clone https://github.com/Olascoaga/aria-omics
cd aria-omics
bash install.sh
```

The installer configures API keys for Anthropic and/or Google (Gemini),
downloads the PBMC 3k test dataset, and verifies the full pipeline.

**Requirements:** Ubuntu / WSL2, Python 3.11, conda or miniforge

---

## Quick Start

```bash
conda activate aria-env
aria
```

```
  ARIA — Agentic Research Intelligence for -omics Analysis

  Action [new/exit]: new
  Data path: /data/h9_experiment
  Your question: I have bulk RNA-seq from human H9 cells.
                 Three conditions: BMAL1 knockout, REV-ERBa knockout,
                 and wildtype, 3 replicates each. What genes are
                 differentially expressed in each knockout?

  CHECKPOINT 1 — Data Audit
  Found: bulk RNA-seq count matrix, 9 samples, Homo sapiens

  CHECKPOINT 2.1 — Experimental Groups
  Detected groups: BMAL1 (3), REV_ERBa (3), WT (3) — confirm?

  ...design walkthrough (organism, factor, batch, replication)...

  CHECKPOINT 2 — Analysis Plan
    Step 1: [bulk_rna_agent] DESeq2 differential expression
    DE thresholds: padj < 0.05, |log2FC| > 0.58
    Contrasts: BMAL1 vs WT, REV_ERBa vs WT, BMAL1 vs REV_ERBa
  Proceed?

  Quality audit passed — no issues found.

  Running analysis...
  → Report: ~/.aria/reports/aria_20260508_.../report.html
```

---

## Adversarial review of every claim

On the validated RNA path, every associative-or-stronger claim is challenged by
a deterministic, LLM-free *devil's advocate* that enumerates the standard
technical confounders (batch effect, ambient RNA, doublets, composition shift,
low replication) and marks — from the run's own structured evidence — which were
addressed and which remain open. Unaddressed alternatives are attached to the
report as caveats. This makes adversarial review reproducible: the same evidence
always yields the same challenges.

For interpretation-heavy, lower-confidence contexts ARIA can additionally run a
multi-round LLM **DebateCouncil** (Proposer vs Critic):

```
Proposer:  "Cluster 3 represents terminally exhausted CD8+ T cells
            (PDCD1+, TOX+, p < 0.001)"

Critic:    "ALTERNATIVE HYPOTHESIS: Could be precursor-exhausted (Tpex).
            TCF7 must be explicitly NEGATIVE, not merely non-significant.
            Request: TCF7 log2FC value."

Verdict:   ACCEPT_REVISED

Consensus: "Cluster 3 shows markers consistent with exhausted CD8+
            T cells (PDCD1+, TOX+). Terminal vs precursor-exhausted
            status requires explicit TCF7 quantification."
```

---

## Development Mode (zero API cost)

```bash
# In ~/.aria/.env
ARIA_DEV_MODE=true
ARIA_DEV_PROVIDER=gemini   # free tier — or "ollama" for local GPU
```

---

## Conda environments

| Environment | Key tools |
|-------------|-----------|
| `aria-rna-env` | scanpy, pydeseq2, gseapy, blitzgsea, scrublet |
| `aria-chromatin-env` | pysam, pybedtools, MACS3, episcanpy, muon |
| `aria-hic-env` | cooler, cooltools, hic-straw, pairtools, chromosight |
| `aria-integration-env` | MOFA+, scGLUE, SCENIC+, decoupler, muon |

For publication or archival runs, generate explicit Linux lockfiles before
tagging or sharing a report:

```bash
scripts/generate_locks.sh
```

Reports embed `envs/*.linux-64.lock` when present and show a visible warning
when lockfiles are missing.

---

## Roadmap

```
v4.0     done     DesignAgent, all pairwise contrasts, paper-theme reports
v4.1     done     AuditAgent — quality linter gates dispatch before DESeq2
v4.2/4.3 done     GEO/SRA connector — analyze public datasets by accession
v4.3.1   done     Hardening: thread-safe bus/memory, prompt cache, exp logs
v4.3.2   done     scRNA end-to-end validated; rna_qc handles raw 10x matrices
v4.3.3   done     scRNA multi-sample workflow: per-sample QC + concat + Harmony
v4.3.4   done     Pseudobulk DE between conditions (rna_pseudobulk_de) +
                  lognorm-counts auto-recovery for Seurat-derived h5ads
v4.3.5   done     NarrativeAgent E2E for scRNA / pseudobulk — UMAP +
                  per-cell-type DE table + pathway dotplots in HTML report
v4.3.6   done     Trajectory analysis (PAGA + root-gated DPT) E2E +
                  adaptive connectivity reporting for mature populations
v4.3.7   done     Cell-cell communication (LIANA) E2E + autocrine
                  exclusion + adaptive rank-metric fallback
v4.3.8   done     TUI/Orchestrator → NarrativeAgent shape normalisation:
                  the path is now end-to-end ready for scRNA reports
v4.3.9   done     Pseudobulk DE between conditions auto-dispatched on
                  comparison intent — TUI now runs aging / treatment
                  contrasts per cell type in one shot
v4.3.10  done     TUI scRNA report hardening — normalized NarrativeAgent
                  inputs and fixed latent report-generation bugs
v4.3.11  done     Production hardening — bulk RNA regression fixed, pytest
                  wrapper added, mocks require explicit dev opt-in
v4.3.12  done     Stability closeout — h5ad obs design inference,
                  processed-h5ad QC, grounded scRNA narrative, scRNA TSV
                  supplements, large-dataset resume/cache guards. Reuses
                  pre-existing obs cell-type annotation when present
                  (skips Leiden/CellTypist), hard-stops on unrecoverable
                  scaled matrices, and surfaces per-cluster DE timeouts as
                  visible warnings on atlas-scale inputs.
v4.4     done     Publication readiness — composition correction, global
                  FDR, ORA backgrounds, power statements, provenance,
                  lockfiles, methodology.json, reproducible mode
v4.5     done     Raw ingestion bridge — deterministic 10X matrix-triplet
                  ingestion and gated FASTQ/kb planning/execution
v4.5.x   done     Reliability, governance & reproducibility hardening before
                  scATAC — deterministic narrative + devil's advocate,
                  planned-vs-run ledger, single-source provenance stamping,
                  design honesty (covariates / explicit contrasts / no
                  filename fallback / propagated thresholds), anti-fabrication
                  guard, complete IPC contracts, blocking 3-tier CI,
                  air-gapped mode
v4.6     done     scATAC release line out of pre-release. scATAC matrix
                  workflow is beta + requires explicit acknowledgement:
                  measured QC, TF-IDF/LSI clustering, per-cluster DA,
                  replicate-gated pseudobulk DA, local motif enrichment, and
                  chromatin narrative blocks. ChIP/CUT&RUN/CUT&TAG remain
                  scaffolded.
post-v4.6 done    scRNA production hardening: raw-count pseudobulk handoff,
                  descriptive marker claims, data-intrinsic QC, root-gated
                  DPT, CellTypist confidence/fallback disclosure, integration
                  overcorrection escalation, 10X MEX directory grouping, and
                  honest QC-failure report blocks.
v4.7     current  Bulk ATAC: QC + MACS3 peak-calling beta slice is open;
                  remaining work is DA via DESeq2 on peak counts
v4.8              IntegrationAgent (WNN + MOFA+ + peak2gene) — deferred
                  until both modalities work standalone
v4.9              Interactive HTML report (sortable tables, plotly figures)
v5.0              Docker image, HPC support, bioRxiv preprint
```

---

## License

MIT — free for academic and commercial use.

---

*Built for biologists who have better things to do than debug pipelines.*
