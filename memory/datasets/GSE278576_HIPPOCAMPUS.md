---
name: GSE278576 hippocampus dataset — Samael's reference
description: Detalles del dataset multi-ómico que Samael usa para validar pipelines de ARIA
type: project
status: active
source_of_truth_for: gse278576_hippocampus_dataset
last_updated: 2026-05-14
originSessionId: 73e36a6b-c405-47b2-ba0f-b480b922999b
---
## El dataset

GSE278576 — Human hippocampus, lifespan study. **Multi-modal: RNA + ATAC** (mismas células / mismos donadores).

Aquí Samael lo tiene preprocesado en Seurat. ARIA hasta ahora solo ha tocado la modalidad RNA:

- **RNA h5ad consolidado**: `/home/medusa/Samael/Single/Samael_Final_Pipeline/data_inputs/reference/GSE278576_hippocampus_RNA.h5ad`
  - 295K cells × ~31K genes
  - 40 donadores (orig.ident)
  - 18 cell-type subclasses (Astro, CA1, CA2-CA3, Chandelier, DG, Endo, LAMP5, Macro, Microglia, NR2F2, OPC, Oligo, PVALB, SUB, SST, T-Cell, VIP, VLMC)
  - obs: age_group (20-39, 40-59, 60-79, 80-100), Gender, percent.mt, percent.ribo, nCount_RNA, nFeature_RNA, subclass, leiden, seurat_clusters
  - obsm: X_pca, X_integrated.rpca, X_umap.rpca
  - raw.X es **log-normalizado** (Seurat), no counts. Auto-recovery en rna_pseudobulk_de via expm1×nCount_RNA/10000.

- **ATAC**: existe pero vive en un **archivo separado** (verificado 2026-05-11: el h5ad RNA NO tiene capa ATAC ni en `layers` ni en `obsm`; solo X_pca / X_umap.rpca / X_integrated.rpca, todos RNA-derived). Samael confirmó **10x Multiome same-cell** — obs_names idénticos cuando se exponga la ruta. v4.4 (scATAC E2E) activará esta modalidad por separado; integración WNN+MOFA+ se mueve a v4.6.

**Ruta del scATAC confirmada (2026-05-11 noche):** Samael tiene un pipeline `Erosion` paralelo en `/home/medusa/Samael/Erosion/` con:

- **40 `.h5mu` paired same-cell** en `data_inputs/muon_processed/<sample>_paired.h5mu` — formato muon multi-omic, RNA + ATAC en el mismo archivo, mismos obs_names. Output del pipeline Erosion `01_fetch_and_pair.py` que descarga GSE278576 de GEO.
- **Subset estratificado** en `results/01_stratified/` (focused on Oligo per `03_qc.py` thresholds: MAX_MITO=5%, MIN_GENES=500).
- Scripts maduros adicionales que ARIA puede aprender / reusar: `04a_motif_scan_oligos.py`, `04b_celloracle_calibracion.py`, `04c_celloracle_inference.py`, `05_cicero_coaccessibility.R`, `05a_GRN_comparison.py`, `05c_structural_analysis.py`, `05d_target_coverage.py`, `05d_topological_architecture.py`.

**Estructura confirmada del h5mu (vía h5py inspection):**

| Modalidad | hc11 shape | Features |
|-----------|------------|----------|
| RNA  | 3,143 × 36,601 | Full transcriptome (no HVG), sparse CSR, gene symbols |
| ATAC | 3,143 × 60,990 | Peaks `chr<n>:start-end`, sparse CSR, Cell Ranger ARC output |

Same-cell pairing confirmado: 3,143 obs_names idénticos entre RNA/ATAC. NO tiene LSI/UMAP precomputado (raw counts del 10x ARC).

**Cell-type labels:** NO vienen en el h5mu; vienen del RNA reference (`GSE278576_hippocampus_RNA.h5ad`, columna `subclass`) por intersection de obs_names. El script Erosion `02_stratify_and_export.py` hace este transfer.

**Decisión validación v4.4 (Samael 2026-05-11):** Usar `hc11_paired.h5mu` (single-donor) replicando el patrón v4.3.2 — validar fast antes de multi-donor.

## Subsets validados en disco (no perder)

- `/tmp/aria_e2e_hippo_multi/out/` — 3 donadores (hc11, hc77, hc1153) procesados de cero con ARIA scRNA pipeline (QC + Harmony + Leiden + CellTypist + DE + ORA). `annotated.h5ad` (455 MB) tiene `X_umap`, `leiden`, `cell_type_celltypist`, `batch`, `sample_id`.

- `/tmp/aria_e2e_hippo_pseudobulk/` — pseudobulk DE del consolidado 40-donor (`pseudobulk_de.csv`, `pathways_per_cluster.csv`, `e2e_report.json`). 11/18 subclases con DE significativo en age_group 80-100 vs 20-39.

- `/tmp/aria_v435_e2e/report/` — primer HTML scRNA publication-ready (5.9 MB) generado por NarrativeAgent en v4.3.5.

## Targets de linaje (para trajectory v4.3.6)

El hipocampo adulto tiene **oligodendrogénesis continua** — el linaje OPC → Oligo es válido para PAGA + DPT incluso en cerebro adulto. Filtrar el h5ad consolidado a:
- `subclass ∈ {OPC, Oligo}` (y opcionalmente Astro como outgroup)
- O directamente del subset multi-sample `/tmp/aria_e2e_hippo_multi/out/annotated.h5ad` donde ya hay clusters anotados

## Para cuando se aborde la modalidad ATAC (v4.4)

El dataset es ideal para validar IntegrationAgent porque RNA y ATAC vienen de las **mismas células** (no diferentes donadores). WNN + peak2gene + MOFA+ pueden correr sobre los mismos `obs_names`. Cuando llegue v4.4: preguntar a Samael por la ruta del h5ad ATAC equivalente.

**Why:** Es el dataset de referencia para todas las validaciones scRNA de ARIA. Saber su estructura ahorra re-exploración.
**How to apply:** Cuando Samael diga "el hipocampo" o "el dataset", asumir este; cuando se necesite un linaje biológicamente válido, usar OPC→Oligo.
