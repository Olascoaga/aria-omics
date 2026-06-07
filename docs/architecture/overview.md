# Architecture Overview

ARIA is not just a collection of pipelines. The central product is a supervised
semantic layer around omics analysis:

1. understand the biological question;
2. inspect the input data;
3. confirm experimental design;
4. run deterministic modality-specific code;
5. record decisions and warnings;
6. generate a report grounded in real output files.

## Main Agent Flow

```mermaid
flowchart TD
    U[User question + local data] --> O[OrchestratorAgent]
    O --> DA[DataAuditAgent]
    DA --> CP1[Checkpoint 1: detected data]
    CP1 --> D[DesignAgent]
    D --> CP2[Design checkpoints: groups, organism, factor, batch, replicates]
    CP2 --> PLAN[Checkpoint 2: analysis plan]
    PLAN -->|confirm| A[AuditAgent]
    PLAN -->|modify thresholds| CP3[Checkpoint 3: thresholds]
    CP3 --> A
    A -->|blocking issues| CP35[Checkpoint 3.5: proceed or cancel]
    A -->|no blocking issues| SETUP[SetupAgent environment check]
    CP35 -->|proceed| SETUP
    SETUP --> DISPATCH[Modality dispatch]
    DISPATCH --> B[BulkRNAAgent]
    DISPATCH --> S[scRNAAgent]
    DISPATCH --> C[ChromatinAgent scaffolded]
    DISPATCH --> H[GenomeArchAgent scaffolded]
    B --> MAYBEINT{2+ modalities or integration requested?}
    S --> MAYBEINT
    C --> MAYBEINT
    H --> MAYBEINT
    MAYBEINT -->|yes| INT[IntegrationAgent scaffolded]
    MAYBEINT -->|no| N[NarrativeAgent]
    INT --> N
    N --> CP5[Checkpoint 5: final review]
    CP5 --> R[HTML report + TSV supplements + methods]
```

The same diagram is stored as [aria_overview.mmd](../diagrams/aria_overview.mmd).
For implementation-impact checks, use the deeper
[Code Dependency Graph](code_graph.md), which maps runtime dependencies,
ownership boundaries, and test anchors. For repository-wide navigation, use the
generated [Graphify map](graphify/README.md): it includes a queryable
`graph.json`, an interactive `graph.html`, and a file-oriented `GRAPH_TREE.html`.
The frozen v4.5 RNA/evidence-governance benchmarking protocol lives in
[benchmarking_v45.md](benchmarking_v45.md); it closes v4.5 as a benchmarking
specification and leaves full benchmark execution as the RNA preprint lane.

## Components

| Component | Responsibility |
|---|---|
| OrchestratorAgent | Owns checkpoint flow and dispatch |
| DataAuditAgent | Detects modalities and data structure |
| DesignAgent | Confirms experimental design before compute |
| AuditAgent | Runs pre-dispatch quality checks |
| SetupAgent | Checks computational environment before modality agents run |
| BulkRNAAgent | Orchestrates bulk RNA count/FASTQ workflows |
| scRNAAgent | Orchestrates QC, integration, clustering, annotation, DE, pseudobulk, beta trajectory and LIANA |
| IntegrationAgent | Runs only for multimodal analyses or when explicitly requested; still scaffolded |
| NarrativeAgent | Writes reports from structured outputs and warnings |
| EnvironmentManager | Runs modality scripts in isolated Conda stacks using JSON IPC |
| ParameterAdvisor | Scores parameter candidates and records decisions |
| ARIAMemory | Stores decisions/findings in SQLite |
| MessageBus | Moves internal agent messages and checkpoint events |

## Script Boundary

Analytical scripts in `aria/scripts/` are subprocess entry points. They should:

- accept a JSON-serializable parameter dict;
- return a JSON-serializable result dict;
- avoid direct access to the message bus;
- validate inputs and return structured errors;
- use explicit mock gates only for development/test modes;
- write output files and summaries that can be validated on resume.

## Dependency Isolation

ARIA separates analytical stacks because bioinformatics libraries often have
conflicting compiled dependencies.

| Stack | Typical tools |
|---|---|
| RNA | Scanpy, AnnData, pyDESeq2, gseapy, LIANA, CellTypist |
| Chromatin | MACS3, pysam, pybedtools, muon / episcanpy |
| Hi-C | cooler, cooltools, hic-straw, pairtools |
| Integration | MOFA+, muon, scGLUE / related integration tools |
