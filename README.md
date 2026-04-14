# ARIA
### Agentic Research Intelligence for -omics Analysis

> *You ask the biological question. ARIA does the rest.*

---

## What is ARIA?

ARIA is an open-source agentic system that automates multi-omics analysis.
Instead of running dozens of tools manually, you describe what you want to
know — and a coordinated team of AI agents handles the rest.

**Supported modalities:**
- scRNA-seq & bulk RNA-seq
- Spatial transcriptomics (Visium, Xenium, MERFISH, Slide-seq)
- scATAC-seq & bulk ATAC-seq
- HiC / Micro-C (3D genome organization)
- ChIP-seq
- CUT&RUN / CUT&TAG
- Multimodal integration of any combination above

**ARIA produces:**
- Quality control reports with biologically-informed thresholds
- Differential analysis with objective parameter justification
- Cross-modal integration with explicit confidence scores
- Narrative reports in plain language
- Fully reproducible methods sections (exportable for manuscripts)

---

## Design Principles

**Local-first** — Your data never leaves your machine.

**Language as interface** — Ask a biological question. ARIA translates it
into an analysis plan, executes it, and explains what it found.

**Supervised autonomy** — Five checkpoints let you review and correct
before critical decisions are made. At Checkpoint 3, you type the
approved parameter value — not just click approve.

**Honest uncertainty** — Every finding carries a confidence level
(HIGH / MEDIUM / LOW / INSUFFICIENT). ARIA tells you when data is
ambiguous or underpowered.

**Institutional memory** — Every approved parameter decision is stored
in a local SQLite database. Over time, ARIA learns your lab's analytical
preferences and cites historical decisions in its justifications.

**Token-efficient** — Inter-agent communication uses CavemanMode
compression. Only user-facing outputs are in normal prose.

**Provider-agnostic** — Works with Claude (Anthropic), Gemini (Google),
or local models via Ollama. Switch providers in one line of config.

---

## Architecture

```
ARIA
├── OrchestratorAgent        Parses questions, designs plans, coordinates
├── DataAuditAgent           Auto-detects data types [CHECKPOINT 1]
├── RNAAgent                 bulk RNA-seq + scRNA-seq + spatial
├── ChromatinAgent           ATAC + ChIP + CUT&RUN + CUT&TAG
├── GenomeArchAgent          HiC, TADs, loops, compartments
├── IntegrationAgent         Multimodal synthesis
└── NarrativeAgent           Reports and visualizations

LLMProvider     Universal LLM abstraction (Anthropic/Gemini/Ollama)
ContextManager  Dynamic context window management per model profile
ParameterAdvisor 3-layer hyperparameter decisions with institutional memory
ARIAMemory      Hierarchical SQLite (Wings/Halls/Rooms/Tunnels)
MessageBus      Inter-agent communication with CavemanMode compression
TUI             Terminal interface (Rich)
```

---

## Installation

```bash
git clone https://github.com/aria-omics/aria
cd aria
bash install.sh
```

The installer guides you through API key configuration for Anthropic
and/or Google (Gemini), downloads the PBMC 3k test dataset, and
verifies the full pipeline.

---

## Quick Start

```bash
conda activate aria-env
aria
```

```
  ARIA -- Agentic Research Intelligence for -omics Analysis

  Action [new/exit]: new
  Data path: /data/lupus_experiment
  Your question: What transcription factors are differentially
                 active in lupus T cells vs healthy controls?

  CHECKPOINT 1 -- Data Audit Results
  ARIA found the following data:
    [+] scRNA-seq: 8 files
    [+] scATAC-seq: 8 files
    Organism: Homo sapiens
    Genome:   hg38

  Is this correct?
  [1] Confirm and continue  [2] Correct metadata  [3] Cancel
```

---

## Development Mode (zero API cost)

```bash
# In ~/.aria/.env
ARIA_DEV_MODE=true
ARIA_DEV_PROVIDER=gemini   # free tier, 1M tokens/day
```

Iterate and debug with Gemini Flash at no cost.
Switch to Claude Sonnet for final production runs.

---

## Public Data Integration

```
Your question: Analyze GSE189903 and compare with my local ATAC data
```

Supported sources: GEO, SRA, ArrayExpress, ENCODE, 4D Nucleome

---

## Checkpoints

| # | When | Question |
|---|------|---------|
| 1 | After data scan | "This is what I found — is it correct?" |
| 2 | Before analysis | "Here is my plan — shall I proceed?" |
| 3 | Parameter decisions | "Type the approved value to confirm" |
| 4 | Mid-analysis | "Here are preliminary findings — adjust?" |
| 5 | Final | "Analysis complete — review report?" |

---

## Citation

If ARIA contributes to a publication, please cite:

```
ARIA: Agentic Research Intelligence for -omics Analysis
[Authors] (2025). bioRxiv. https://doi.org/[TBD]
```

---

## License

MIT — free for academic and commercial use.

---

*Built for biologists who have better things to do than debug pipelines.*
