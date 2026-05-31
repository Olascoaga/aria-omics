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

    DISPATCH --> ICP[Internal parameter checkpoints]
    ICP --> BASE[aria/agents/base_agent.py publish_blocking_escalation]
    BASE --> BUS[aria/bus/message_bus.py wait_for_checkpoint_resolution]
    BUS --> TUI_CP[aria/tui.py or aria/headless.py resolves checkpoint]
    TUI_CP --> ORCH_CP[orchestrator_agent.on_checkpoint_resolved]
    ORCH_CP --> BUS
    BUS --> ICP

    DISPATCH --> RAW[aria/agents/raw_ingestion_agent.py]
    RAW --> RAW_UTILS[aria/utils/raw_ingestion.py]
    RAW_UTILS --> CANONICAL[canonical input files and provenance]
    CANONICAL --> DISPATCH

    DISPATCH --> BULK_AGENT[aria/agents/bulk_rna_agent.py]
    BULK_AGENT --> ENV[aria/utils/environment_manager.py]
    ENV --> BULK_SCRIPT[aria/scripts/rna_bulk_de.py]
    BULK_SCRIPT --> CLASSIFIER[aria/utils/count_classifier.py raw-count guard]
    BULK_SCRIPT --> PATHWAY_VIZ[aria/scripts/rna_pathway_viz.py]
    BULK_SCRIPT --> BULK_OUT[bulk findings: QC, DE, ORA, GSEA, figures, tables]

    DISPATCH --> SCRNA_AGENT[aria/agents/scrna_agent.py]
    SCRNA_AGENT --> BASE
    SCRNA_AGENT --> SCRNA_QC[aria/scripts/rna_qc.py]
    SCRNA_AGENT --> SCRNA_CLUSTER[aria/scripts/rna_clustering.py]
    SCRNA_AGENT --> SCRNA_DA[aria/scripts/rna_diff_abundance.py quasi-Poisson DA]
    SCRNA_AGENT --> SCRNA_PB[aria/scripts/rna_pseudobulk_de.py]
    STATS[aria/utils/stats.py shared statistical helpers] --> SCRNA_DA
    STATS --> SCRNA_PB
    SCRNA_PB --> CLASSIFIER
    SCRNA_AGENT --> SCRNA_PW[aria/scripts/rna_pathway_per_cluster.py]
    SCRNA_PB -->|per-cluster ORA universe| SCRNA_PW
    SCRNA_AGENT --> SCRNA_CCC[aria/scripts/rna_cellcomm.py]
    SCRNA_AGENT --> SCRNA_TRAJ[aria/scripts/rna_trajectory.py]
    SCRNA_QC --> SCRNA_OUT[scRNA findings]
    SCRNA_CLUSTER --> SCRNA_OUT
    SCRNA_DA --> SCRNA_OUT
    SCRNA_PB --> SCRNA_OUT
    SCRNA_PW --> SCRNA_OUT
    SCRNA_CCC --> SCRNA_OUT
    SCRNA_TRAJ --> SCRNA_OUT

    AUDIT --> MUDATA[aria/utils/mudata_io.py .h5mu real reader]
    DISPATCH --> CHROM[aria/agents/chromatin_agent.py scaffolded]
    CHROM --> CHROM_QC[aria/scripts/chromatin_qc.py measured-only QC]
    CHROM_QC --> MUDATA
    DISPATCH --> HIC[aria/agents/genome_arch_agent.py scaffolded]
    DISPATCH --> INT_GATE[Integration validation gate]
    INT_GATE --> INT[aria/agents/integration_agent.py scaffolded]
    INT --> BASE

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
    NARR --> LLM[aria/llm/provider.py]
    LLM --> LLM_USAGE[llm_usage.jsonl deterministic controls and cost]
    LLM_USAGE --> METHODOLOGY
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
| `rna_pseudobulk_de.py` or `rna_diff_abundance.py` | `scrna_agent.py`, `ScrnaNarrator`, FDR-strategy and power report text | These scripts drive publication-facing inferential claims. DA uses quasi-Poisson (overdispersion-corrected, ADR-018); do not revert to plain Poisson — it anti-conservatively gates the composition covariate. |
| `rna_pseudobulk_de.py` per-group `background_genes` / `scrna_agent` `background_genes_by_cluster` / `rna_pathway_per_cluster.py` universe | per-cluster ORA significance, `ScrnaNarrator` Methods + per-block background size | The ORA universe is per cell type (genes tested in that cluster's pseudobulk, ADR-018). A single global background inflates per-cluster enrichment; legacy results without per-group background fall back to the global universe. |
| `rna_pseudobulk_de.py` composition covariate / `composition_skipped_reason` | `ScrnaNarrator` composition caveats, `scrna_agent` composition gating | The self-proportion composition covariate is dropped when collinear with the contrast (`_abs_corr` >= `COMPOSITION_COLLINEARITY_MAX`, C3/ADR-021); do not re-enable it unconditionally — it inflates variance exactly when abundance shifts with condition. |
| `bus/message_bus.py` indices / persistence / `get_pending_checkpoints(experiment_id=...)` | TUI/headless polls, `OrchestratorAgent.run` (`enable_persistence`), checkpoint resolution, crash recovery | Findings/escalations are served from eviction-consistent indices and optionally persisted per-run (R6/ADR-021). Keep `_index`/`_deindex` in sync with the FIFO deque, keep persistence per-`experiment_id`, and scope reads by experiment so concurrent runs sharing the global bus don't cross-read. |
| `stats.py` BH correction helper | `rna_pseudobulk_de.py`, `rna_diff_abundance.py`, FDR tests | Shared multiple-testing code must stay numerically stable; duplicating local BH implementations risks divergent significance calls. |
| `count_classifier.py` raw-count detection | `rna_bulk_de.py` (`_load_counts` hard-refuse), `rna_pseudobulk_de.py` (integer-likeness + log-norm recovery probes), `count_source` provenance | The single detector deciding raw vs normalized input for DESeq2. Loosening `is_raw_counts` lets normalized matrices become pseudo-counts; the sampler is seeded — keep it deterministic for reproducible mode. |
| `message_bus.py`, `BaseAgent.publish_blocking_escalation`, or checkpoint handling in `tui.py` / `headless.py` | `scrna_agent.py` Leiden resolution, `integration_agent.py` WNN/MOFA, `orchestrator_agent.py` CP3 handling | Internal parameter checkpoints must block script execution until user/headless resolution; otherwise custom/skip choices are decorative. |
| `orchestrator_agent.py` CP3 resolution | CP3 threshold tuning, internal agent parameter checkpoints, dispatch thread lifecycle | Internal CP3 messages carry `agent_parameter_checkpoint=True` and must not trigger threshold-tuning redispatch. |
| `orchestrator_agent.py` integration validation gate | `IntegrationAgent`, multimodal report sections, registry integrity | Scaffolded WNN/MOFA+/peak-to-gene code must not dispatch until validation is closed; otherwise beta scripts can emit publication-looking integration output. |
| `parameter_advisor.py` metric evaluators | scRNA clustering CP3, WNN k CP3, memory decisions | Candidate metrics shown to users must be measured or explicitly marked as not computed; fabricated WNN weights and zero-filled modularity are invalid. |
| `scrna_agent.py` focus or annotation fallback logic | DesignIntelligence, focused h5ad materialization, report labels | Runtime logic must match explicit obs values or external annotation output; do not reintroduce tissue/cell-type alias maps or hardcoded marker panels under ADR-011. |
| `NarrativeBlock` schema | every modality narrator, validators, renderer, `methodology.json` | This is the report evidence contract. |
| `validators.py` | all report generation | Validators are the last integrity gate before claims reach HTML. The causal guard scans ARIA's authored claim, not external named entities (DB term names, gene symbols) carried in evidence; `collect_named_entities` is also reused by the render-level prose scan. |
| `compose_prose.py` or `render_blocks.py` | HTML findings for all block-backed modalities | Rendering changes can turn valid results into cryptic or misleading reports. |
| `llm/provider.py` or `utils/provenance.py` LLM usage schema | `NarrativeAgent` report provenance, `methodology.json`, prompt cache behavior | Narrative confidence/prose must remain reproducible: deterministic controls, model tier, token counts, and cache semantics are part of audit provenance. Every call is time-bounded (`timeout`, R3) and tier fallbacks are recorded as degradation (`is_fallback`/`fallback_*` → `collect_llm_usage` `degraded`/`fallback_calls`, R4/ADR-020) — keep these fields flowing to the report. |
| `environment_manager.py` | all script-running agents | It controls conda stack execution and JSON IPC boundaries. Scripts run under `Popen(start_new_session=True)`; a timeout reaps the whole process group via `_terminate_process_tree`/`os.killpg` (R5/ADR-020) — do not revert to bare `subprocess.run`, which orphans BLAS/numba grandchildren. `_resolve_env` is preferred-env-if-installed → FALLBACK_ENV with NO aliasing (B12); do not reintroduce an `env_aliases.json` read — no writer exists and it contradicts SetupAgent's "no aliases" policy. |
| `chromatin_qc.py` metric helpers / `mudata_io.py` `.h5mu` reader | `chromatin_agent.py`, v4.6 scATAC QC, DataAudit `.h5mu` detection | Chromatin QC must emit only measured metrics (ADR-019/ADR-002): TSS/FRiP/barcodes are real or `None`+`metrics_not_computed`, never placeholders. `.h5mu` is the detected scATAC entry; the MuData reader returns structured blockers when tooling/ATAC modality is absent — do not fabricate. |
| `data_audit_agent.py` SIGNATURES / `_scan_directory` extensions | modality classification, CP1, `.h5mu`/chromatin routing | `.h5mu` (paired RNA+ATAC) is scanned and classified as `scATAC` (C8); changing the signature order or extension set can make the v4.6 entry input undetectable again. |
| module-global accessors `bus` / `env_manager` | all agents and tests importing global coordination helpers | These are lazy accessors, not eagerly constructed singletons; avoid import-time workspace creation or broker state just by importing enums/classes. |
| `memory.py` decisions schema | checkpoints, provenance, reports, resume | DB shape changes can break old sessions and audit trails. Operational `memory/` files are private local context and must stay ignored/out of GitHub. |

## Ownership Boundaries

Agents own orchestration and user-visible decisions. Scripts own deterministic
analysis. Narrators own translation from structured results to
`NarrativeBlock`. Renderers own presentation.

Do not cross these boundaries casually:

- Scripts should not publish bus messages or write report prose.
- Agents should not parse generated HTML to recover results.
- Agents that publish in-dispatch parameter checkpoints must use
  `publish_blocking_escalation`; fire-and-forget checkpoints are invalid for
  user-controllable parameters.
- Narrators should not invent analyses that are absent from structured output.
- scRNA focus/annotation code must use explicit obs values, external annotation
  output, or unresolved labels; no built-in tissue marker panels or dataset
  aliases.
- Renderers should not create scientific claims that are not present in a
  block.
- Validators should downgrade or block unsafe claims, not silently rewrite
  biological meaning.

## Impact Checklist

Before changing a production path, answer these questions:

1. What result keys does this module produce or consume?
2. Which checkpoint text or confirmed design fields will change?
3. If it publishes a checkpoint from the dispatch thread, what blocks execution
   until `on_checkpoint_resolved` records a user/headless decision?
4. Which report sections, figures, tables, or `methodology.json` fields depend
   on it?
5. Does the change affect replay/resume or old workspaces?
6. Does it touch publication claims: DE, global/local FDR, power, ORA/GSEA,
   LIANA, PAGA/DPT, provenance, or dependency locks?
7. Which focused tests prove the contract, and which smoke test covers the
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
- Raw-count guard (DESeq2 input integrity): `tests/test_count_classifier.py`,
  `tests/test_bulk_raw_count_guard.py`
- Chromatin v4.6 readiness (no fabricated QC, `.h5mu` detection/reader, no env
  alias path): `tests/test_chromatin_readiness.py` (mudata-gated cases skip
  without the chromatin stack)
- LLM/subprocess reliability (call timeout, model-degradation provenance,
  current model IDs, process-group kill on timeout):
  `tests/test_llm_reliability.py`
- Bus durability + per-run isolation + indexed reads:
  `tests/test_bus_durability.py`
- Composition-covariate collinearity guard:
  `tests/test_composition_collinearity.py` (end-to-end case is pydeseq2-gated)
- B7/B11 housekeeping and ADR-011 guards:
  `tests/test_pytest_smoke.py::test_global_bus_and_env_manager_are_lazy_accessors`,
  `tests/test_pytest_smoke.py::test_marker_fallback_annotation_is_explicit_and_conservative`,
  `tests/test_pytest_smoke.py::test_scrna_infers_and_materializes_cell_focus`,
  `tests/test_pytest_smoke.py::test_scrna_cell_focus_does_not_expand_without_explicit_focus`,
  `tests/test_pytest_smoke.py::test_design_intelligence_scrna_focused_group_feasibility`
- Checkpoint blocking and dispatch safety:
  `tests/test_pytest_smoke.py::test_internal_parameter_checkpoint_blocks_until_user_resolution`,
  `tests/test_pytest_smoke.py::test_wnn_checkpoint_skip_prevents_script_execution`,
  `tests/test_pytest_smoke.py::test_orchestrator_does_not_dispatch_on_internal_cp3_resolution`
- Integration scaffold gate and parameter honesty:
  `tests/test_pytest_smoke.py::test_orchestrator_skips_scaffolded_integration_agent`,
  `tests/test_pytest_smoke.py::test_wnn_advice_does_not_fabricate_pre_run_metrics`,
  `tests/test_pytest_smoke.py::test_wnn_checkpoint_marks_pre_run_metrics_as_not_computed`,
  `tests/test_pytest_smoke.py::test_leiden_subprocess_modularity_is_not_replaced_with_zero`,
  `tests/test_registry_integrity.py::test_scaffold_integration_agent_is_not_dispatched`
- LLM deterministic provenance:
  `tests/test_pytest_smoke.py::test_llm_provider_forces_deterministic_generation`,
  `tests/test_pytest_smoke.py::test_llm_cache_key_includes_deterministic_controls`,
  `tests/test_pytest_smoke.py::test_collect_llm_usage_summarizes_deterministic_provenance`,
  `tests/test_pytest_smoke.py::test_provenance_section_renders_llm_usage`
- Main integration smoke: `tests/test_pytest_smoke.py`
