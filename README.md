# ARIA

### Agentic Research Intelligence for -omics Analysis

> *You ask the biological question. ARIA does the rest.*

![Version](https://img.shields.io/badge/version-4.1-blue)
![Status](https://img.shields.io/badge/status-active-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.11-blue)

---

## What is ARIA?

ARIA is an open-source agentic system that automates multi-omics analysis
for biological researchers. Instead of running dozens of tools manually,
you describe what you want to know — and a coordinated team of AI agents
handles the rest.

**Supported modalities:**

- scRNA-seq and bulk RNA-seq
- scATAC-seq and bulk ATAC-seq
- ChIP-seq (TF binding and histone marks)
- CUT&RUN / CUT&TAG
- HiC / Micro-C (3D genome: TADs, loops, compartments A/B)
- Multimodal integration of any combination above

**ARIA produces:**

- Pre-analysis quality audit with actionable warnings before expensive compute
- Differential expression with all pairwise contrasts (not just vs control)
- Pathway enrichment (ORA + GSEA) per contrast
- Publication-ready figures: volcano plots, PCA/MDS, ORA dotplots, GSEA running sums
- HTML report in paper/publication style with embedded figures
- Manuscript-ready methods section with exact parameters used
- Reproducible decision log (every threshold choice stored in SQLite)

---

## Current Status — May 2026

```
INFRASTRUCTURE                     STATUS
────────────────────────────────────────────────────
MessageBus + CavemanMode           done
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
BulkRNAAgent (DESeq2)              done  ✓ H9 BMAL1/REV-ERBα
  rna_bulk_de.py                   done
  rna_pathway_viz.py               done
scRNAAgent                         done  ✓ PBMC 3k
  rna_qc.py  (MAD + stress ctx)    done
  rna_clustering.py (Leiden)       done
ChromatinAgent                     done
  chromatin_qc.py                  done
  chromatin_peaks.py (MACS3)       done
GenomeArchAgent                    done
  hic_inspect.py                   done
  hic_qc_and_balance.py            done
  hic_topology.py (out-of-core)    done
NarrativeAgent (HTML report)       done  ✓ paper theme
IntegrationAgent (WNN + MOFA+)     building
GEO/SRA connectors                 planned
```

**End-to-end validated** on human H9 cells (bulk RNA-seq, 3 conditions × 3 replicates):

- All pairwise contrasts generated: BMAL1 vs WT, REV-ERBα vs WT, BMAL1 vs REV-ERBα
- DESeq2 + ORA (GO_BP, KEGG, Reactome) + GSEA per contrast
- Publication-ready HTML report with embedded volcano plots and pathway dotplots

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

**Institutional memory** — Every approved parameter decision is stored
in a local SQLite database. Over time, ARIA learns your lab's analytical
preferences and cites historical decisions in its justifications.

**Dependency isolation** — Each analytical stack runs in its own Conda
environment. IPC via JSON files prevents C-library conflicts between
scanpy, cooler, MACS3, and pysam.

**Token-efficient** — Inter-agent communication uses CavemanMode
compression. Only user-facing outputs are in normal prose.

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
  AuditAgent              Quality linter before dispatch   [CHECKPOINT 3.5 if blocking]
    · Replicate correlation (outlier / swap detection)
    · PCA batch dominance
    · STAR alignment rate
  BulkRNAAgent            DESeq2, all pairwise contrasts, ORA, GSEA
  scRNAAgent              QC, clustering, annotation, DE
  ChromatinAgent          ATAC + ChIP + CUT&RUN + CUT&TAG
  GenomeArchAgent         HiC, TADs, loops, compartments A/B
  IntegrationAgent        Multimodal synthesis (WNN, MOFA+)  [building]
  NarrativeAgent          HTML report + methods section

  LLMProvider       Universal LLM abstraction (Anthropic / Gemini / Ollama)
  ContextManager    4-step degradation cascade for local models
  ParameterAdvisor  3-layer hyperparameter decisions + institutional memory
  EnvironmentManager IPC via JSON, isolated Conda stacks per modality
  DebateCouncil     Internal peer review: Proposer vs Critic (2–3 rounds)
  ARIAMemory        Hierarchical SQLite: Wings / Halls / Rooms / Findings
  MessageBus        Inter-agent pub/sub with CavemanMode compression
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

---

## Roadmap

```
v4.0  done     DesignAgent, all pairwise contrasts, paper-theme reports
v4.1  done     AuditAgent — quality linter gates dispatch before DESeq2
v4.2  next     GEO/SRA connector — analyze public datasets by accession
v4.3           IntegrationAgent (WNN + MOFA+), multimodal reports
v4.4           Interactive HTML report (sortable tables, plotly figures)
v5.0           Docker image, HPC support, bioRxiv preprint
```

---

## License

MIT — free for academic and commercial use.

---

*Built for biologists who have better things to do than debug pipelines.*
