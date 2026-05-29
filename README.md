# ARIA

### Agentic Research Intelligence for -omics Analysis

> *You ask the biological question. ARIA does the rest.*

![Version](https://img.shields.io/badge/version-4.5.4-blue)
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

- **Production-like validated:** bulk RNA-seq count matrices; scRNA-seq
  single-sample and multi-sample workflows; processed `.h5ad` pseudobulk
  workflows with design metadata in `obs`.
- **Validated / beta:** bulk RNA-seq FASTQ preprocessing, trajectory
  summaries, LIANA cell-cell communication, processed `.h5ad` recovery, and
  GEO/SRA connector paths.
- **Scaffolded / roadmap:** scATAC-seq, bulk ATAC-seq, ChIP-seq, CUT&RUN /
  CUT&TAG, full Hi-C / Micro-C workflows, and multimodal WNN/MOFA+
  integration.

**ARIA produces:**

- Pre-analysis quality audit with actionable warnings before expensive compute
- Differential expression with all pairwise contrasts (not just vs control)
- Pathway enrichment (ORA + GSEA) per contrast
- Publication-ready figures: volcano plots, PCA/MDS, ORA dotplots, GSEA running sums
- HTML report in paper/publication style with embedded figures
- Manuscript-ready methods section with exact parameters used
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

## Current Status — May 2026

```
INFRASTRUCTURE                     STATUS
────────────────────────────────────────────────────
MessageBus + compact wire format   done
ARIAMemory (SQLite hierarchical)   done
LLMProvider + ContextManager       done
ParameterAdvisor (3-layer)         done
EnvironmentManager (IPC/JSON)      done
DebateCouncil (peer review)        done

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
  rna_trajectory.py (PAGA+DPT)     beta — validated on hippocampus subset
  rna_cellcomm.py (LIANA)          beta — validated on GSE278576
ChromatinAgent                     scaffolded
  chromatin_qc.py                  done
  chromatin_peaks.py (MACS3)       scaffolded
GenomeArchAgent                    scaffolded
  hic_inspect.py                   done
  hic_qc_and_balance.py            done
  hic_topology.py (out-of-core)    scaffolded
NarrativeAgent (HTML report)       done  ✓ paper theme
  Narrative kernel (evidence cards) done — validated NarrativeBlock objects
  Claim Compiler (evidence tiers)   done — per-claim tier + manifest (X14)
IntegrationAgent (WNN + MOFA+)     scaffolded — pending end-to-end validation
GEO/SRA connectors                 done   ✓ GSE183948 validated
```

The current stable baseline is `v4.5.4`: bulk RNA + scRNA core paths are
closed for practical use, publication-readiness provenance is embedded in
reports, raw-ingestion planning/conversion is available for supported 10X
inputs, and reports are composed from validated modality blocks for scRNA and
bulk RNA.

`v4.5.3` tagged the pre-ATAC integrity freeze: centralized version metadata,
installer secret hygiene, registry-integrity checks, scaffold dispatch gating,
typed script IPC contracts, design-matrix validation before DESeq2, synthetic
ground-truth DE benchmarking, scientific QC red flags, and the **Claim
Compiler**. `v4.5.4` adds scientific-honesty hardening: pseudobulk scRNA now
defaults to per-cluster FDR while still reporting global FDR, power is reported
against the effective global-BH threshold, and log-normalized count recovery is
visibly low-confidence. Pseudobulk significant-gene counts can differ from
pre-`v4.5.4` reports by design.

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
  tables, LIANA interactions, and PAGA/DPT trajectory summaries. The legacy
  cell-level `rna_de_per_cluster.py` stage may time out on atlas-scale inputs;
  this is non-blocking when pseudobulk outputs are available.

**Publication-grade HTML report for scRNA / pseudobulk** (v4.3.5) — run
`tests/test_scrna_e2e.py --emit-html` and ARIA writes a self-contained
HTML to `<workspace>/report/` containing: UMAP coloured by cell-type and
condition, per-cell-type DE summary bar, per-(group×comparison) DE table
sorted by significance with top up/down genes, ORA dotplots for the top
cell types across GO_BP / KEGG / Reactome, Methods section with exact
thresholds and design formula, and the parameter-decisions log.

**Trajectory analysis (PAGA + DPT)** (v4.3.6) — preprocessed h5ads can
be dispatched with `--trajectory-h5ad PATH --trajectory-groupby COL
--trajectory-root TYPE`. Validated on a 50k-cell OPC + Oligo + Astro
subset from the hippocampus dataset: DPT pseudotime correctly orders
OPC (0.216) → Oligo (0.239) → Astro (0.367) with OPC selected as root.
PAGA reports `max_connectivity = 0.00231` — characteristic of mature /
non-developmental tissue — and the HTML embeds both a normal and a
log-scaled cluster graph so weak adult-tissue edges remain visible
without the report making false developmental claims. RNA velocity is
automatically skipped (the dataset lacks spliced / unspliced layers);
re-quantification via velocyto or kb-python `nac` mode is documented
in the Methods section as the path to enable it.

**Cell-cell communication (LIANA)** (v4.3.7) — dispatch with
`--cellcomm-h5ad PATH --cellcomm-groupby cell_type_celltypist`.
Validated on the GSE278576 multi-sample annotated h5ad (9 cell types):
LIANA `rank_aggregate` with n_perms=50 runs in 7 s. Autocrine pairs
(source == target) are excluded a priori (~1.5k of ~13.5k rows on this
dataset). The script auto-falls back to `specificity_rank` when LIANA
emits all-NaN `magnitude_rank` (a recent-version quirk). Top non-
autocrine hits recovered are classical glia-neuron signaling axes:
APOE → TREM2 (Astro → Micro), C3 → NRP1 (Micro → L2-3 neurons),
VCAN → TLR2 (OPC → Micro), NRG1 → MS4A4A (InN VIP → Micro). HTML
includes a sender × receiver interaction-count heatmap, the top-N
L-R bar chart, and a sortable table with CellPhone p-values.

---

## Design Principles

**Local-first** — Your data never leaves your machine.

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
  ChromatinAgent          ATAC + ChIP + CUT&RUN + CUT&TAG  [scaffolded]
  GenomeArchAgent         HiC, TADs, loops, compartments   [scaffolded]
  IntegrationAgent        Conditional multimodal synthesis  [scaffolded]
  NarrativeAgent          HTML report + methods section
  Final Review            Accept / revise / export          [CHECKPOINT 5]

  LLMProvider       Universal LLM abstraction (Anthropic / Gemini / Ollama)
  ContextManager    4-step degradation cascade for local models
  ParameterAdvisor  3-layer hyperparameter decisions + institutional memory
  EnvironmentManager IPC via JSON, isolated Conda stacks per modality
  DebateCouncil     Internal peer review: Proposer vs Critic (2–3 rounds)
  ARIAMemory        Hierarchical SQLite: Wings / Halls / Rooms / Findings
  MessageBus        Inter-agent pub/sub with compact internal messages
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

## The DebateCouncil in action

Every biological interpretation with MEDIUM or LOW confidence goes
through internal peer review before reaching the user:

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
v4.3.6   done     Trajectory analysis (PAGA + DPT) E2E + adaptive
                  connectivity reporting for mature populations
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
v4.6     next     scATAC end-to-end — chromatin_agent + chromatin_qc +
                  chromatin_peaks already scaffolded; need LSI clustering
                  + differential accessibility + motifs + narrative module
v4.7              ATAC bulk end-to-end — DA via DESeq2 on peak counts
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
