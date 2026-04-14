# ARIA 🧬

### Agentic Research Intelligence for -omics Analysis

> \*You ask the biological question. ARIA does the rest.\*

!\[Tests](https://img.shields.io/badge/tests-96%2F96%20passing-brightgreen)
!\[Status](https://img.shields.io/badge/status-alpha%20v0.1-blue)
!\[License](https://img.shields.io/badge/license-MIT-green)
!\[Python](https://img.shields.io/badge/python-3.11-blue)

\---

## What is ARIA?

ARIA is an open-source agentic system that automates multi-omics analysis
for biological researchers. Instead of running dozens of tools manually,
you describe what you want to know — and a coordinated team of AI agents
handles the rest.

**Supported modalities:**

* scRNA-seq and bulk RNA-seq (including spatial transcriptomics)
* scATAC-seq and bulk ATAC-seq
* ChIP-seq (TF binding and histone marks)
* CUT\&RUN / CUT\&TAG
* HiC / Micro-C (3D genome: TADs, loops, compartments A/B)
* Multimodal integration of any combination above

**ARIA produces:**

* Quality control reports with biologically-informed thresholds
* Differential analysis with objective parameter justification
* 3D genome topology (compartments, TADs, chromatin loops)
* Narrative reports in plain language
* Fully reproducible methods sections — exportable for manuscripts

\---

## Current Status — April 2026

```
INFRASTRUCTURE                     STATUS    TESTS
──────────────────────────────────────────────────
MessageBus + CavemanMode           done      21/21
ARIAMemory (SQLite hierarchical)   done
LLMProvider + ContextManager       done
ParameterAdvisor (3-layer)         done
EnvironmentManager (IPC/JSON)      done      15/15
DebateCouncil (peer review)        done      15/15

AGENTS                             STATUS    TESTS
──────────────────────────────────────────────────
DataAuditAgent                     done
OrchestratorAgent                  done
RNAAgent (bulk + scRNA)            done
  rna\_qc.py  (MAD + stress ctx)    done      live PBMC 3k
  rna\_clustering.py (Leiden)       done
ChromatinAgent                     done      24/24
  chromatin\_qc.py                  done
  chromatin\_peaks.py (MACS3)       done
GenomeArchAgent                    done      21/21
  hic\_inspect.py                   done
  hic\_qc\_and\_balance.py            done
  hic\_topology.py (out-of-core)    done

IN PROGRESS
──────────────────────────────────────────────────
IntegrationAgent (WNN + MOFA+)     building
NarrativeAgent (HTML/PDF report)   planned
GEO/SRA connectors                 planned
Spatial transcriptomics module     planned

TOTAL: 96/96 tests green
```

**End-to-end validated** on PBMC 3k (10x Genomics, 2,700 human PBMCs):

* QC: 2,700 to 2,643 cells, adaptive MAD thresholds
* Clustering: ParameterAdvisor recommended resolution=0.4, silhouette=0.720
* Annotation: Claude Haiku identified T cells, B cells, NK cells, Monocytes
* DebateCouncil example: revised "terminally exhausted CD8+ T cells" claim
to require explicit TCF7 quantification (Tex vs Tpex distinction)

\---

## Design Principles

**Local-first** — Your data never leaves your machine.

**Language as interface** — Ask a biological question. ARIA translates it
into an analysis plan, executes it, and explains what it found.

**Supervised autonomy** — Five checkpoints let you review and correct
before critical decisions are made. At Checkpoint 3, you type the
approved parameter value — active confirmation, not rubber-stamping.

**Honest uncertainty** — Every finding carries a confidence level
(HIGH / MEDIUM / LOW / INSUFFICIENT). ARIA tells you when data is
ambiguous or underpowered.

**Institutional memory** — Every approved parameter decision is stored
in a local SQLite database. Over time, ARIA learns your lab's analytical
preferences and cites historical decisions in its justifications.

**Dependency isolation** — Each analytical stack runs in its own Conda
environment (aria-rna-env, aria-chromatin-env, aria-hic-env,
aria-integration-env). IPC via JSON files prevents C-library conflicts
between scanpy, cooler, MACS3, and pysam.

**Token-efficient** — Inter-agent communication uses CavemanMode
compression. Only user-facing outputs are in normal prose.

**Provider-agnostic** — Works with Claude (Anthropic), Gemini (Google),
or local models via Ollama. Switch providers in one config line.

\---

## Architecture

```
ARIA
  OrchestratorAgent       Parses questions, designs plans, coordinates
  DataAuditAgent          Auto-detects data types  \[CHECKPOINT 1]
  RNAAgent                bulk RNA-seq + scRNA-seq + spatial
  ChromatinAgent          ATAC + ChIP + CUT\&RUN + CUT\&TAG
  GenomeArchAgent         HiC, TADs, loops, compartments A/B
  IntegrationAgent        Multimodal synthesis (WNN, MOFA+)  \[in progress]
  NarrativeAgent          Reports and visualizations          \[planned]

  LLMProvider       Universal LLM abstraction (Anthropic/Gemini/Ollama)
  ContextManager    4-step degradation cascade for local models
  ParameterAdvisor  3-layer hyperparameter decisions + institutional memory
  EnvironmentManager IPC via JSON, isolated Conda stacks per modality
  DebateCouncil     Internal peer review: Proposer vs Critic (2-3 rounds)
  ARIAMemory        Hierarchical SQLite: Wings/Halls/Rooms/Findings/Tunnels
  MessageBus        Inter-agent communication with CavemanMode compression
  TUI               Terminal interface (Rich)
```

\---

## Installation

```bash
git clone https://github.com/Olascoaga/aria-omics
cd aria-omics
bash install.sh
```

The installer configures API keys for Anthropic and/or Google (Gemini),
downloads the PBMC 3k test dataset, and verifies the full pipeline.

**Requirements:** Ubuntu / WSL2, Python 3.11, conda or miniforge

\---

## Quick Start

```bash
conda activate aria-env
aria
```

```
  ARIA -- Agentic Research Intelligence for -omics Analysis

  Action \[new/exit]: new
  Data path: /data/lupus\_experiment
  Your question: What transcription factors are differentially
                 active in lupus T cells vs healthy controls?

  CHECKPOINT 1 -- Data Audit Results
  ARIA found the following data:
    \[+] scRNA-seq: 8 files
    \[+] scATAC-seq: 8 files
    Organism: Homo sapiens
    Genome:   hg38

  Is this correct?
  \[1] Confirm and continue  \[2] Correct metadata  \[3] Cancel
```

\---

## The DebateCouncil in action

Every biological interpretation with MEDIUM or LOW confidence goes
through internal peer review before reaching the user:

```
Proposer:  "Cluster 3 represents terminally exhausted CD8+ T cells
            (PDCD1+, TOX+, p < 0.001)"

Critic:    "ALTERNATIVE HYPOTHESIS: Could be precursor-exhausted (Tpex),
            not terminal Tex. TCF7 must be explicitly NEGATIVE,
            not merely non-significant.
            Request: TCF7 log2FC value."

Verdict:   ACCEPT\_REVISED

Consensus: "Cluster 3 shows markers consistent with exhausted CD8+
            T cells (PDCD1+, TOX+). Terminal vs precursor-exhausted
            status requires explicit TCF7 quantification."

Limitations:
  - TCF7 status requires quantification for Tex/Tpex distinction
  - Small cluster size (n=87) limits statistical power
  - Protein-level validation recommended before publishing
```

This example demonstrates how ARIA prevents overclaiming in the
manuscript — a distinction that matters in immunology and that
automated pipelines typically miss.

\---

## Development Mode (zero API cost)

```bash
# In \~/.aria/.env
ARIA\_DEV\_MODE=true
ARIA\_DEV\_PROVIDER=gemini   # free tier — or "ollama" for local GPU
```

Develop and debug with Gemini Flash at no cost.
Switch to Claude Sonnet for final production runs ($0.01-0.05 per analysis).

\---

## Conda environments

Each analytical stack is isolated to prevent C-library conflicts:

|Environment|Key tools|
|-|-|
|`aria-rna-env`|scanpy, pydeseq2, squidpy, scrublet, gseapy|
|`aria-chromatin-env`|pysam, pybedtools, MACS3, episcanpy, muon|
|`aria-hic-env`|cooler, cooltools, hic-straw, pairtools, chromosight|
|`aria-integration-env`|MOFA+, scGLUE, SCENIC+, decoupler, muon|

```bash
conda env create -f envs/aria-rna-env.yml
conda env create -f envs/aria-chromatin-env.yml
conda env create -f envs/aria-hic-env.yml
conda env create -f envs/aria-integration-env.yml
```

\---

## Checkpoints

|#|When|What ARIA asks|
|-|-|-|
|1|After data scan|"This is what I found — is it correct?"|
|2|Before analysis|"Here is my plan — shall I proceed?"|
|3|Parameter decisions|Type the approved value to confirm|
|4|Mid-analysis|"Here are preliminary findings — adjust?"|
|5|Final|"Analysis complete — review report?"|

Checkpoint 3 requires typing the approved value (not clicking a button)
to ensure active engagement with critical analytical decisions.

\---

## Run the validation tests

```bash
conda activate aria-env

# Core infrastructure
python tests/test\_integration.py          # 21/21

# EnvironmentManager + scripts
python tests/test\_environment\_manager.py  # 15/15

# DebateCouncil peer review
python tests/test\_debate\_council.py       # 15/15

# ChromatinAgent
python tests/test\_chromatin\_agent.py      # 24/24

# GenomeArchAgent
python tests/test\_genome\_arch\_agent.py    # 21/21

# End-to-end with real PBMC 3k data
python tests/test\_pbmc\_e2e.py
```

\---

## Roadmap

```
v0.1  current   RNA + Chromatin + 3D genome, DebateCouncil, 96 tests
v0.2  next      IntegrationAgent (WNN + MOFA+), NarrativeAgent
v0.3            GEO/SRA connectors, spatial transcriptomics
v0.4            Benchmark vs manual analysis, bioRxiv preprint
v1.0            Production-ready, Docker image, HPC support
```

## Citation

```
ARIA: Agentic Research Intelligence for -omics Analysis
\[Authors] (2025). bioRxiv. https://doi.org/\[TBD]
```

\---

## License

MIT — free for academic and commercial use.

\---

*Built for biologists who have better things to do than debug pipelines.*

