---
name: Paths importantes en la máquina de Samael
description: Ubicaciones de repo, datos, genomas, reportes y configuración — para no preguntar dónde están
type: reference
status: active
source_of_truth_for: local_paths
last_updated: 2026-05-14
originSessionId: 30506a0f-929d-4e63-be74-d496dc9b0be5
---
Máquina: **Medusa** (WSL2 / Ubuntu).

## Repositorio y configuración

| Path | Qué es |
|------|--------|
| `~/Samael/ARIA/` | Repo principal (`github.com/Olascoaga/aria-omics`) |
| `~/.aria/config.yaml` | Config de LLM provider |
| `~/.aria/.env` | Private environment file; never copy or commit its contents. |

## Datos y procesamiento

| Path | Qué es |
|------|--------|
| `~/Samael/H9-RNA/raw_fastq/` | 18 FASTQs reales, 9 samples (B1-B3 / R1-R3 / WT1-WT3) — dataset de validación bulk RNA-seq |
| `~/Samael/H9-RNA/aria_processing/{qc,aligned,counts}/` | Outputs intermedios del pipeline bulk |
| `~/aria-data/pbmc3k_test/{matrix.mtx,barcodes.tsv,genes.tsv}` | pbmc3k MEX — dataset canónico para validación scRNA (descargado por installer) |
| `~/Samael/Single/prueba/GSE278576_hc*_raw_feature_bc_matrix.h5` | 40 muestras de hipocampo humano (GSE278576), raw 10x CellRanger output |
| `~/Samael/Erosion/raw/GSE278576_hippocampus_RNA.h5ad` | h5ad consolidado de hipocampo (Seurat-procesado, 295K cells × 2K HVGs, ya con `subclass` y `orig.ident`) |
| `~/Samael/Single/Samael_Final_Pipeline/data_inputs/reference/GSE278576_hippocampus_RNA.h5ad` | **El h5ad de referencia** del usuario para análisis de aging. Mismo contenido que arriba pero esta es la copia "oficial" en su pipeline final. 2.2 GB. `raw.X` está log-normalizado (Seurat NormalizeData scale_factor=10000), `pseudobulk_de` lo revierte automáticamente con `nCount_RNA` |

## Recursos compartidos

| Path | Qué es |
|------|--------|
| `~/.aria/genomes/hg38/{genome.fa, annotation.gtf, star_index/}` | Genoma hg38 + GTF + STAR index (preparado por SetupAgent) |
| `~/.aria/reports/aria_<timestamp>_<slug>_<id>/` | Cada reporte vive en su propia carpeta con `report.html`, `figures/`, `tables/` |
| `~/.aria/memory.db` | SQLite jerárquica (Wings/Halls/Rooms/Findings/Decisions/Tunnels) |
| `~/.aria/llm_cache/` | Prompt cache (v4.3.1+, deshabilitable con `ARIA_LLM_CACHE=0`) |
| `~/.aria/logs/exp_<id>.log` | Log persistente por experimento (v4.3.1+) |
| `~/.aria/workspace/` | Inputs/outputs JSON de scripts en ejecución |
| `~/.aria/workspace/failed/<stack>_<run_id>/` | Postmortem de runs fallidos (v4.3.1+) |

## Dataset de validación bulk

H9 BMAL1 KO / REV-ERBα KO (9 samples) — el experimento canónico end-to-end:
- 94.3% unique mapping
- 1489 shared DE genes (Jaccard 0.4 entre BMAL1 KO y REV-ERBα KO)
- Top genes: IGFBP5, COL3A1, CDKN1A, ACAT2
- Hallazgo: convergencia en p53 signaling entre ambos KOs (potencial paper)
