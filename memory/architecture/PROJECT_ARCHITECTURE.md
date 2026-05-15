---
name: ARIA — Arquitectura y mapa de archivos clave
description: Dónde vive cada componente, qué hace, cuántas líneas tiene
type: project
status: active
source_of_truth_for: repository_architecture
last_updated: 2026-05-14
originSessionId: 992731fd-4099-4636-a948-1a891d380475
---
## Directorio raíz: /home/medusa/Samael/ARIA

```
aria/
  tui.py                          (818)  — TUI Rich, polling loop de checkpoints en background
  __init__.py                     (4)    — package version

  agents/
    base_agent.py                 (170)  — clase base abstracta con run()
    orchestrator_agent.py         (677)  — parsea pregunta → plan → coordina agentes
    data_audit_agent.py           (418)  — auto-detecta tipos de dato [CHECKPOINT 1]
    design_agent.py               (566)  — diseño experimental interactivo [CP 2.1–2.6]
    audit_agent.py                (436)  — quality linter antes de dispatch [CP 3.5]
    bulk_rna_agent.py             (507)  — DESeq2, ORA, GSEA
    scrna_agent.py                (827)  — QC, clustering, anotación, DE single-cell
    chromatin_agent.py            (637)  — ATAC, ChIP, CUT&RUN, CUT&TAG
    genome_arch_agent.py          (756)  — HiC, TADs, loops, compartimentos A/B
    integration_agent.py          (741)  — WNN + MOFA+ (scaffoldeado, pendiente validación)
    narrative_agent.py            (1947) — genera HTML report paper-style [el más grande]
    debate_council.py             (409)  — peer review interno: Proposer vs Critic
    setup_agent.py                (719)  — configuración inicial de entorno

  scripts/                        — pipelines bioinformáticos (corren en envs conda aislados)
    rna_bulk_de.py                (2258) — DESeq2 pipeline completo [el más grande]
    rna_quantify.py               (554)
    rna_fastq_qc.py               (357)
    rna_align.py                  (329)
    rna_qc.py                     (182)
    rna_clustering.py             (90)
    rna_pathway_viz.py            (378)
    rna_integration.py            (105)  — multi-sample RNA
    rna_trajectory.py             (165)
    rna_cellcomm.py               (178)
    chromatin_qc.py               (417)
    chromatin_peaks.py            (295)  — MACS3
    hic_inspect.py                (143)
    hic_qc_and_balance.py         (304)
    hic_topology.py               (511)  — out-of-core
    integration_wnn.py            (294)  — Weighted Nearest Neighbor
    integration_mofa.py           (330)  — MOFA+
    integration_peak2gene.py      (448)  — peak-to-gene links
    _base.py                      (152)  — base para scripts

  connectors/
    geo_connector.py              (515)  — descarga GEO (GSE/SRP/PRJNA), infiere diseño

  llm/
    provider.py                   (318)  — abstracción Anthropic/Gemini/Ollama
    context_manager.py            (402)  — cascada de degradación de 4 pasos
    parameter_advisor.py          (654)  — 3 capas de decisión + memoria institucional

  bus/
    message_bus.py                (129)  — pub/sub inter-agente con CavemanMode

  memory/
    memory.py                     (270)  — SQLite jerárquico (Wings/Halls/Rooms/Findings)

  utils/
    environment_manager.py        (347)  — IPC via JSON, stacks conda aislados

tests/                            (~4076 líneas totales)
  test_bulk_rna.py
  test_scrna.py
  test_chromatin_agent.py
  test_genome_arch_agent.py
  test_integration.py / test_integration_agent.py
  test_narrative_agent.py
  test_debate_council.py
  test_environment_manager.py
  test_pbmc_e2e.py
```

## Entornos Conda

| Entorno | Herramientas clave |
|---------|-------------------|
| `aria-rna-env` | scanpy, pydeseq2, gseapy, blitzgsea, scrublet |
| `aria-chromatin-env` | pysam, pybedtools, MACS3, episcanpy, muon |
| `aria-hic-env` | cooler, cooltools, hic-straw, pairtools, chromosight |
| `aria-integration-env` | MOFA+, scGLUE, SCENIC+, decoupler, muon |

## Flujo principal (checkpoints)

```
DataAuditAgent [CP1] → DesignAgent [CP 2.1–2.6] → Plan [CP2]
  → AuditAgent [CP3.5 si bloqueante]
  → Dispatcher (background thread)
    → BulkRNA / scRNA / Chromatin / GenomeArch / Integration
  → NarrativeAgent → HTML report [CP5]
```

**Why:** Mapa de referencia para navegar el código sin tener que re-explorar.
**How to apply:** Usar como índice cuando se necesita editar un componente específico.
