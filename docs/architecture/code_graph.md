# Code Dependency Graph

This document is the working impact map for ARIA code changes. It is not a
complete Python import graph. It tracks the runtime data flow, ownership
boundaries, and production-impact edges that matter when adding features or
polishing existing behavior.

Use it before changing shared agents, script result schemas, report rendering,
or checkpoint/design logic.

## RNA And Reporting Flow

```mermaid
flowchart TD
    TUI[aria/tui.py] --> GEO[aria/connectors/geo_connector.py]
    TUI --> ORCH[aria/agents/orchestrator_agent.py]
    GEO --> AUDIT[aria/agents/data_audit_agent.py]
    ORCH --> AUDIT

    AUDIT --> CP1[Checkpoint 1 detected data]
    CP1 --> DESIGN[aria/agents/design_agent.py]
    DESIGN --> CP2[Checkpoint 2 design and plan]
    CP2 --> DI[aria/agents/design_intelligence.py]
    DI --> AUDIT2[aria/agents/audit_agent.py]
    AUDIT2 --> SETUP[aria/agents/setup_agent.py]
    SETUP --> DISPATCH[Modality dispatch]

    DISPATCH --> RAW[aria/agents/raw_ingestion_agent.py]
    RAW --> RAW_UTILS[aria/utils/raw_ingestion.py]
    RAW_UTILS --> CANONICAL[canonical input files and provenance]
    CANONICAL --> DISPATCH

    DISPATCH --> BULK_AGENT[aria/agents/bulk_rna_agent.py]
    BULK_AGENT --> ENV[aria/utils/environment_manager.py]
    ENV --> BULK_SCRIPT[aria/scripts/rna_bulk_de.py]
    BULK_SCRIPT --> PATHWAY_VIZ[aria/scripts/rna_pathway_viz.py]
    BULK_SCRIPT --> BULK_OUT[bulk findings: QC, DE, ORA, GSEA, figures, tables]

    DISPATCH --> SCRNA_AGENT[aria/agents/scrna_agent.py]
    SCRNA_AGENT --> SCRNA_QC[aria/scripts/rna_qc.py]
    SCRNA_AGENT --> SCRNA_CLUSTER[aria/scripts/rna_clustering.py]
    SCRNA_AGENT --> SCRNA_DA[aria/scripts/rna_diff_abundance.py]
    SCRNA_AGENT --> SCRNA_PB[aria/scripts/rna_pseudobulk_de.py]
    SCRNA_AGENT --> SCRNA_PW[aria/scripts/rna_pathway_per_cluster.py]
    SCRNA_AGENT --> SCRNA_CCC[aria/scripts/rna_cellcomm.py]
    SCRNA_AGENT --> SCRNA_TRAJ[aria/scripts/rna_trajectory.py]
    SCRNA_QC --> SCRNA_OUT[scRNA findings]
    SCRNA_CLUSTER --> SCRNA_OUT
    SCRNA_DA --> SCRNA_OUT
    SCRNA_PB --> SCRNA_OUT
    SCRNA_PW --> SCRNA_OUT
    SCRNA_CCC --> SCRNA_OUT
    SCRNA_TRAJ --> SCRNA_OUT

    DISPATCH --> CHROM[aria/agents/chromatin_agent.py scaffolded]
    DISPATCH --> HIC[aria/agents/genome_arch_agent.py scaffolded]
    DISPATCH --> INT[aria/agents/integration_agent.py scaffolded]

    BULK_OUT --> NARR[aria/agents/narrative_agent.py]
    SCRNA_OUT --> NARR
    CHROM --> NARR
    HIC --> NARR
    INT --> NARR

    NARR --> REGISTRY[aria/agents/narrative/registry.py]
    REGISTRY --> BULK_NARR[aria/agents/narrative/narrators/bulk_rna.py]
    REGISTRY --> SCRNA_NARR[aria/agents/narrative/narrators/scrna.py]
    SCRNA_NARR --> SCRNA_HELPERS[aria/agents/_narrative_scrna.py]
    BULK_NARR --> BLOCKS[NarrativeBlock list]
    SCRNA_NARR --> BLOCKS
    BLOCKS --> VALIDATORS[aria/agents/narrative/validators.py]
    VALIDATORS --> PROSE[aria/agents/narrative/compose_prose.py]
    PROSE --> RENDER[aria/agents/narrative/render_blocks.py]
    RENDER --> REPORT[report.html]
    NARR --> METHODOLOGY[methodology.json]
```

## High-Impact Edges

| If you change | Check these dependents | Why it is risky |
|---|---|---|
| `geo_connector.py` inferred design schema | `data_audit_agent.py`, `design_agent.py`, `bulk_rna_agent.py`, TUI GEO flow | GEO groups, organism, aliases, and sample IDs seed the entire design path. |
| `data_audit_agent.py` h5ad or GEO design inference | `design_agent.py`, `scrna_agent.py`, `design_intelligence.py`, CP1 text | A wrong condition/replicate/groupby guess can run valid code on the wrong biological design. |
| `design_agent.py` design dict | `bulk_rna_agent.py`, `scrna_agent.py`, `design_intelligence.py`, methods/provenance | `groups`, `main_factor`, covariates, aliases, and pseudobulk settings are consumed downstream. |
| `bulk_rna_agent.py` | `rna_bulk_de.py`, narrative bulk narrator, methods, release validation | It maps confirmed design onto count matrix columns and owns bulk result shape. |
| `rna_bulk_de.py` result schema | `BulkRnaNarrator`, `NarrativeAgent`, bulk workflow docs, tests | Report claims, figures, tables, ORA, GSEA, and methodology depend on specific keys. |
| `rna_pathway_viz.py` | `rna_bulk_de.py`, bulk narrative, report artifacts | It creates GSEA/ORA figures and tables that reports reference. |
| `scrna_agent.py` result schema | `_narrative_scrna.py`, `ScrnaNarrator`, scRNA workflow docs | scRNA reports assume stable keys for QC, composition, pseudobulk, pathways, LIANA, and trajectory. |
| `rna_pseudobulk_de.py` or `rna_diff_abundance.py` | `scrna_agent.py`, `ScrnaNarrator`, FDR-strategy and power report text | These scripts drive publication-facing inferential claims. |
| `NarrativeBlock` schema | every modality narrator, validators, renderer, `methodology.json` | This is the report evidence contract. |
| `validators.py` | all report generation | Validators are the last integrity gate before claims reach HTML. |
| `compose_prose.py` or `render_blocks.py` | HTML findings for all block-backed modalities | Rendering changes can turn valid results into cryptic or misleading reports. |
| `environment_manager.py` | all script-running agents | It controls conda stack execution and JSON IPC boundaries. |
| `memory.py` decisions schema | checkpoints, provenance, reports, resume | DB shape changes can break old sessions and audit trails. |

## Ownership Boundaries

Agents own orchestration and user-visible decisions. Scripts own deterministic
analysis. Narrators own translation from structured results to
`NarrativeBlock`. Renderers own presentation.

Do not cross these boundaries casually:

- Scripts should not publish bus messages or write report prose.
- Agents should not parse generated HTML to recover results.
- Narrators should not invent analyses that are absent from structured output.
- Renderers should not create scientific claims that are not present in a
  block.
- Validators should downgrade or block unsafe claims, not silently rewrite
  biological meaning.

## Impact Checklist

Before changing a production path, answer these questions:

1. What result keys does this module produce or consume?
2. Which checkpoint text or confirmed design fields will change?
3. Which report sections, figures, tables, or `methodology.json` fields depend
   on it?
4. Does the change affect replay/resume or old workspaces?
5. Does it touch publication claims: DE, global/local FDR, power, ORA/GSEA,
   LIANA, PAGA/DPT, provenance, or dependency locks?
6. Which focused tests prove the contract, and which smoke test covers the
   integration path?

## Current Narrative Rule

For modalities backed by `NarrativeBlock`, prose is the primary report
presentation. Structured evidence tables, figures, and file links are audit
support. A report that only exposes claim rows and evidence tables is not
considered sufficiently narrative.

## Test Anchors

Use these tests as impact anchors when editing the graph's major nodes:

- Narrative kernel: `tests/test_narrative_types.py`,
  `tests/test_narrative_validators.py`, `tests/test_narrative_render_blocks.py`
- Bulk narrator and GSEA surfacing: `tests/test_narrator_bulk.py`,
  `tests/test_pathway_viz.py`
- scRNA narrator: `tests/test_narrator_scrna.py`
- GEO/design mapping: `tests/test_geo_design.py`
- Main integration smoke: `tests/test_pytest_smoke.py`
