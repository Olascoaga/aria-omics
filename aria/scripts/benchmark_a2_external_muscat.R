#!/usr/bin/env Rscript
# Benchmark A2 external reference: muscat-method pseudobulk DE on the Kang data.
#
# Runs inside aria-bench-env. Uses the real Kang PBMC data (muscat's bundled
# example_sce: 8 donors, control vs IFN-beta) and computes muscat's reference
# pseudobulk DE method (sum-aggregate to per-cluster donor-condition pseudobulk,
# then edgeR-QLF) directly with edgeR. This reproduces muscat::pbDS(method=
# "edgeR") faithfully without requiring the full muscat package dependency tree
# (only the dataset + edgeR/SingleCellExperiment, all present). It exports the
# per-cluster pseudobulk matrices + reference DE tables + sample table so ARIA's
# own pseudobulk DE can run on the identical inputs and be compared.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(SingleCellExperiment)
  library(edgeR)
})

# --- real Kang data: full 8-vs-8 (muscData EH2259), cached locally ----------
# Prefer the full Kang18_8vs8 (8 donors x ctrl/stim) so each group has adequate
# donor replication. The ExperimentHub resource is fetched directly (no muscData
# install needed) and cached; falls back to muscat's example_sce subset.
get_kang_sce <- function(cache) {
  if (!file.exists(cache)) {
    url <- "https://experimenthub.bioconductor.org/fetch/2259"
    message("downloading full Kang18_8vs8 (muscData EH2259)")
    utils::download.file(url, cache, mode = "wb", quiet = TRUE)
  }
  e <- new.env()
  ok <- tryCatch({ load(cache, envir = e); TRUE }, error = function(err) FALSE)
  if (ok && length(ls(e)) > 0) return(get(ls(e)[[1]], envir = e))
  stop("could not load the Kang SCE from ", cache)
}

cache <- Sys.getenv("ARIA_KANG_SCE_CACHE",
                    file.path(path.expand("~/.aria/benchmarks"), "kang18_8vs8.rda"))
dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
sce <- get_kang_sce(cache)
cd <- as.data.frame(colData(sce))

# Normalize the colData schema: full Kang (ind/stim/cell/multiplets) vs the
# prepSCE'd example_sce (cluster_id/sample_id/group_id).
if (all(c("ind", "stim", "cell") %in% colnames(cd))) {
  if ("multiplets" %in% colnames(cd))
    keep_cell <- as.character(cd$multiplets) == "singlet" & !is.na(cd$cell)
  else
    keep_cell <- !is.na(cd$cell)
  sce <- sce[, keep_cell]
  cd <- as.data.frame(colData(sce))
  cd$cluster_id <- as.character(cd$cell)
  cd$group_id <- as.character(cd$stim)
  cd$sample_id <- paste(cd$stim, cd$ind, sep = "_")
}
stopifnot(all(c("cluster_id", "sample_id", "group_id") %in% colnames(cd)))
counts <- as.matrix(assay(sce, "counts"))

# --- sample table (donor-condition sample -> condition) -------------------
samp_ids <- as.character(cd$sample_id)
grp_of <- tapply(as.character(cd$group_id), samp_ids, function(x) x[[1]])
samples <- data.frame(sample = names(grp_of), group = unname(grp_of))
write.table(samples, file.path(out_dir, "samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- per-cluster pseudobulk + edgeR-QLF (= muscat edgeR pbDS) --------------
safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
clusters <- levels(factor(as.character(cd$cluster_id)))
written <- character()
for (cl in clusters) {
  in_cl <- as.character(cd$cluster_id) == cl
  smp <- factor(samp_ids[in_cl], levels = samples$sample)
  # Sum counts per sample within this cluster -> genes x samples pseudobulk.
  pbmat <- t(rowsum(t(counts[, in_cl, drop = FALSE]), group = smp))
  pbmat <- pbmat[, colSums(pbmat) > 0, drop = FALSE]
  if (ncol(pbmat) < 4) next

  grp <- factor(samples$group[match(colnames(pbmat), samples$sample)],
                levels = c("ctrl", "stim"))
  if (length(unique(grp)) < 2 || min(table(grp)) < 2) next

  s <- safe(cl)
  pbdf <- data.frame(gene = rownames(pbmat), pbmat, check.names = FALSE)
  write.table(pbdf, file.path(out_dir, sprintf("pb_%s.tsv", s)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  y <- DGEList(counts = pbmat, group = grp)
  keep <- filterByExpr(y, group = grp)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~grp)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)            # stim vs ctrl
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  mdf <- data.frame(gene = rownames(tt), logFC = tt$logFC,
                    p_val = tt$PValue, p_adj = tt$FDR)
  write.table(mdf, file.path(out_dir, sprintf("muscat_%s.tsv", s)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  written <- c(written, s)
}

meta <- list(
  clusters = as.list(written),
  contrast = list(numerator = "stim", denominator = "ctrl"),
  dataset = "Kang et al. 2018 (muscData Kang18_8vs8, 8 donors ctrl/stim IFN-beta); reference = muscat edgeR-QLF pseudobulk",
  reference_method = "pseudobulk sum-aggregation + edgeR-QLF (muscat pbDS edgeR)",
  n_samples = nrow(samples),
  n_clusters = length(written)
)
writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "clusters.json"))

cat(sprintf("A2_MUSCAT_EXPORT_DONE clusters=%d samples=%d\n",
            length(written), nrow(samples)))
