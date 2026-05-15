---
name: v4.6 IntegrationAgent plan — WNN + MOFA+ + peak2gene
description: Plan detallado guardado para cuando se aborde la integración multimodal RNA+ATAC. Confirmado por Samael 2026-05-11 que se difiere hasta tener scATAC (v4.4) y ATAC bulk (v4.5) cerrados.
type: project
status: active
source_of_truth_for: v46_integration_plan
last_updated: 2026-05-14
originSessionId: 73e36a6b-c405-47b2-ba0f-b480b922999b
---
## Por qué se difiere

Samael decidió 2026-05-11 que ARIA debe poder analizar **cada modalidad por separado** antes de intentar integrarlas. Eso pone:
- v4.4 → scATAC E2E
- v4.5 → ATAC bulk E2E
- v4.6 → este plan (integración multimodal)

## Pre-requisitos confirmados

- Dataset: hipocampo GSE278576 — **10x Multiome same-cell** (confirmado por Samael).
- `obs_names` idénticos entre RNA y ATAC → **WNN aplica directo sin Seurat anchors / MultiVI / GLUE**.
- Para v4.6, asumir que `chromatin_agent` ya produjo un h5ad scATAC procesado (peaks anotados, LSI, clusters) y el RNA h5ad ya tiene QC + Harmony + Leiden + CellTypist (lo que v4.3.3 / v4.3.5 ya emiten).

## Scripts ya scaffolded (1813 LOC sin validar)

| Archivo | LOC | Función |
|---------|-----|---------|
| `aria/scripts/integration_wnn.py` | 294 | Weighted Nearest Neighbors — joint kNN ponderado por modality contribution |
| `aria/scripts/integration_mofa.py` | 330 | MOFA+ factors compartidos vs específicos por modality |
| `aria/scripts/integration_peak2gene.py` | 448 | Correlación peak-gene con vecindad espacial |
| `aria/agents/integration_agent.py` | 741 | Agent wrapper — orquesta los 3 scripts |

## Orden de trabajo

1. **Smoke test cada script independiente** (mismo patrón de v4.3.6/7: validar antes de integrar al narrativo).
   - `integration_wnn.py` primero (más simple, joint embedding).
   - `integration_mofa.py` segundo (variance decomposition).
   - `integration_peak2gene.py` tercero (necesita peak coordinates + gene coordinates).
2. **Adapter:** extender `rna_narrative_adapter.py` con `findings.integration` mapping. O crear `narrative_integration.py` si se vuelve grande.
3. **Módulo narrativo:** `_narrative_integration.py` con:
   - Joint UMAP (clusters WNN vs RNA-only vs ATAC-only — side by side).
   - WNN weight density per cluster (heatmap o stacked bar mostrando contribución RNA/ATAC por cluster).
   - MOFA+ variance plot (top factors × variance explained, marcados por modality).
   - Top peak2gene scatter (correlation vs distance, con anotación de top hits).
4. **Harness:** `--integration-rna PATH --integration-atac PATH` (asume same-cell, mismos obs_names).
5. **Validar `IntegrationAgent` ya escrito** vs reusarlo como standalone scripts:
   - Si el agent existing tiene bus / checkpoints útiles → reusar.
   - Si está hardcoded a un flow → mejor scripts standalone.

## Done criteria v4.6

- HTML self-contained con:
  - Joint UMAP comparison panel
  - WNN cluster-level modality weights
  - MOFA+ top factors interpretados (cell cycle factor flag + variance explained per modality)
  - peak2gene top correlaciones positivas + negativas (poised regulatory elements)
  - Methods section con WNN k, MOFA+ n_factors, peak2gene window size
- Validado sobre el par RNA + ATAC del hipocampo de Samael.

## Caveats anticipados

- **MOFA+ es lento** (~30 min sobre 295K cells × 2K genes × 100K peaks). Vamos a querer subsamplear como hicimos con trayectorias.
- **peak2gene window default es 500kb**. Para hipocampo (heterogéneo), considerar 250kb.
- **Cell cycle factor flag:** MOFA+ a menudo produce 1-2 factors que capturan ciclo celular — el módulo narrativo debe detectarlos (correlación con G2M/S scores si están disponibles) y marcarlos como "excluidos de interpretación biológica".

**Why:** Plan completo guardado para no re-pensarlo cuando se aborde.
**How to apply:** Cuando se complete v4.5 y se arranque v4.6, leer este archivo primero.
