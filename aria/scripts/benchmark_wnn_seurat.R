#!/usr/bin/env Rscript
# V48 Slice I2 — Seurat FindMultiModalNeighbors external reference driver.
# Runs in aria-bench-signac-env (Seurat v5). Consumes the IDENTICAL shared
# barcodes, RNA PCA embedding (X_pca.tsv), and depth-filtered ATAC LSI embedding
# (X_lsi.tsv) exported by benchmark_wnn_export.py, executes canonical Seurat WNN
# (`FindMultiModalNeighbors` with `l2.norm = TRUE` + `FindClusters` on `wsnn`),
# and writes cell-level modality weights, clusters, and weighted.nn neighbors.
#
# Honesty contract (ADR-002): any missing R dependency or runtime error is
# emitted as an explicit `not_run` / `error` status JSON, never fabricated.

suppressWarnings(suppressMessages({
  ok <- requireNamespace("optparse", quietly = TRUE)
}))
if (!ok) {
  cat('{"status":"not_run","reason":"optparse unavailable"}\n')
  quit(status = 0)
}
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--in-dir", type = "character"),
  make_option("--out-cells", type = "character"),
  make_option("--out-knn", type = "character"),
  make_option("--out-json", type = "character"),
  make_option("--k", type = "integer", default = 20L),
  make_option("--resolution", type = "double", default = 0.5),
  make_option("--seed", type = "integer", default = 0L)
)))

write_json <- function(obj, path) {
  parts <- vapply(names(obj), function(k) {
    v <- obj[[k]]
    if (is.numeric(v)) sprintf('"%s": %s', k, format(v, scientific = FALSE))
    else sprintf('"%s": "%s"', k, gsub('"', "'", as.character(v)))
  }, character(1))
  writeLines(paste0("{\n  ", paste(parts, collapse = ",\n  "), "\n}"), path)
}

fail <- function(reason) {
  if (!is.null(opt$`out-json`)) {
    write_json(list(status = "not_run", reason = reason), opt$`out-json`)
  }
  cat(sprintf('{"status":"not_run","reason":"%s"}\n', reason))
  quit(status = 0)
}

for (pkg in c("Matrix", "Seurat")) {
  if (!requireNamespace(pkg, quietly = TRUE)) fail(paste0(pkg, " unavailable"))
}

suppressWarnings(suppressMessages({
  library(Matrix)
  library(Seurat)
}))

indir <- opt$`in-dir`
res <- tryCatch({
  set.seed(opt$seed)
  barcodes <- readLines(file.path(indir, "barcodes.tsv"))
  n_cells <- length(barcodes)

  pca <- as.matrix(read.table(file.path(indir, "X_pca.tsv"), sep = "\t", header = FALSE))
  lsi <- as.matrix(read.table(file.path(indir, "X_lsi.tsv"), sep = "\t", header = FALSE))
  rownames(pca) <- barcodes
  colnames(pca) <- paste0("PC_", seq_len(ncol(pca)))
  rownames(lsi) <- barcodes
  colnames(lsi) <- paste0("LSI_", seq_len(ncol(lsi)))

  # Minimal sparse scaffold required by CreateSeuratObject; WNN operates purely
  # on the precomputed pca and lsi DimReduc embeddings.
  dummy_counts <- Matrix::sparseMatrix(
    i = seq_len(min(n_cells, 4L)),
    j = seq_len(min(n_cells, 4L)),
    x = 1,
    dims = c(4L, n_cells),
    dimnames = list(paste0("feat_", 1:4), barcodes)
  )
  obj <- CreateSeuratObject(counts = dummy_counts, assay = "RNA")
  obj[["pca"]] <- CreateDimReducObject(embeddings = pca, key = "PC_", assay = "RNA")
  obj[["lsi"]] <- CreateDimReducObject(embeddings = lsi, key = "LSI_", assay = "RNA")

  k_eff <- min(as.integer(opt$k), max(2L, n_cells - 2L))
  knn_range <- min(200L, max(k_eff + 1L, n_cells - 1L))

  obj <- FindMultiModalNeighbors(
    object = obj,
    reduction.list = list("pca", "lsi"),
    dims.list = list(seq_len(ncol(pca)), seq_len(ncol(lsi))),
    k.nn = k_eff,
    knn.range = knn_range,
    l2.norm = TRUE,
    modality.weight.name = c("RNA.weight", "ATAC.weight"),
    verbose = FALSE
  )
  obj <- FindClusters(
    object = obj,
    graph.name = "wsnn",
    resolution = opt$resolution,
    random.seed = opt$seed,
    verbose = FALSE
  )

  cells_df <- data.frame(
    barcode = barcodes,
    seurat_cluster = as.character(obj@meta.data$seurat_clusters),
    rna_weight = as.numeric(obj@meta.data$RNA.weight),
    atac_weight = as.numeric(obj@meta.data$ATAC.weight),
    stringsAsFactors = FALSE
  )
  write.table(cells_df, opt$`out-cells`, sep = "\t", quote = FALSE, row.names = FALSE)

  nn_mat <- obj@neighbors$weighted.nn@nn.idx
  nb_strings <- vapply(seq_len(n_cells), function(i) {
    idxs <- nn_mat[i, ]
    idxs <- idxs[idxs != i & idxs >= 1L & idxs <= n_cells]
    paste(barcodes[idxs], collapse = ",")
  }, character(1))
  knn_df <- data.frame(
    barcode = barcodes,
    neighbors = nb_strings,
    stringsAsFactors = FALSE
  )
  write.table(knn_df, opt$`out-knn`, sep = "\t", quote = FALSE, row.names = FALSE)

  list(
    status = "success",
    seurat_version = as.character(packageVersion("Seurat")),
    n_cells = n_cells,
    n_clusters = length(unique(cells_df$seurat_cluster)),
    mean_rna_weight = round(mean(cells_df$rna_weight), 4),
    mean_atac_weight = round(mean(cells_df$atac_weight), 4),
    k = k_eff
  )
}, error = function(e) {
  list(status = "error", reason = conditionMessage(e))
})

if (!is.null(opt$`out-json`)) write_json(res, opt$`out-json`)
cat(sprintf('{"status":"%s"}\n', res$status))
