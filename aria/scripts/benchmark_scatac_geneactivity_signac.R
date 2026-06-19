#!/usr/bin/env Rscript
# scATAC P4.3 — Signac GeneActivity reference driver for the gene-activity concordance
# benchmark. Runs in aria-bench-signac-env. Builds the canonical Signac gene-activity
# estimate from the SAME fragments + cells ARIA scored: counts fragments over each
# gene's body + promoter (gene body extended `extend-upstream` bp upstream, the Signac
# GeneActivity convention) via CreateFragmentObject + FeatureMatrix, then writes the
# per-gene MEAN activity (gene, score) that the pure-numpy P4.3 scorer consumes.
#
# No fabrication: any missing dependency / asset is an honest non-success status in the
# output JSON, never an invented gene-activity matrix (ADR-002 discipline).

suppressWarnings(suppressMessages({ ok <- requireNamespace("optparse", quietly = TRUE) }))
if (!ok) { cat('{"status":"not_run","reason":"optparse unavailable"}\n'); quit(status = 0) }
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--fragments", type = "character"),
  make_option("--barcodes", type = "character"),
  make_option("--gtf", type = "character"),
  make_option("--out-scores", type = "character"),
  make_option("--out-json", type = "character"),
  make_option("--extend-upstream", type = "double", default = 2000),
  make_option("--cores", type = "integer", default = 1)
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
  if (!is.null(opt$`out-json`)) write_json(list(status = "not_run", reason = reason), opt$`out-json`)
  cat(sprintf('{"status":"not_run","reason":"%s"}\n', reason)); quit(status = 0)
}

need <- c("Signac", "Matrix", "GenomicRanges", "rtracklayer", "IRanges", "S4Vectors")
for (pkg in need) if (!requireNamespace(pkg, quietly = TRUE)) fail(paste0(pkg, " unavailable"))
suppressWarnings(suppressMessages({
  library(Signac); library(Matrix); library(GenomicRanges); library(rtracklayer)
}))

# Parallelize FeatureMatrix over gene regions (single-threaded is slow on a
# whole-genome fragments file). Signac uses the future framework.
if (!is.null(opt$cores) && opt$cores > 1 && requireNamespace("future", quietly = TRUE)) {
  future::plan("multicore", workers = opt$cores)
  options(future.globals.maxSize = 8 * 1024^3)
}

res <- tryCatch({
  barcodes <- readLines(opt$barcodes)

  # Gene body + promoter ranges from the local GTF (gene_name keyed), extended upstream.
  gtf <- rtracklayer::import(opt$gtf)
  genes <- gtf[gtf$type == "gene"]
  genes <- genes[!is.na(genes$gene_name) & genes$gene_name != ""]
  genes <- genes[!duplicated(genes$gene_name)]
  coords <- Signac::Extend(genes, upstream = opt$`extend-upstream`, downstream = 0)
  names(coords) <- genes$gene_name

  frags <- CreateFragmentObject(path = opt$fragments, cells = barcodes, validate.fragments = FALSE)
  ga <- FeatureMatrix(fragments = frags, features = coords, cells = barcodes)

  # FeatureMatrix keys rows by the range string (chr-start-end), not the gene; map back
  # to gene_name so the symbols match ARIA's gene scores (GeneActivity does the same).
  range_str <- GRangesToString(coords, sep = c("-", "-"))
  gene_of <- setNames(names(coords), range_str)
  mapped <- gene_of[rownames(ga)]
  ga <- ga[!is.na(mapped), , drop = FALSE]
  rownames(ga) <- mapped[!is.na(mapped)]

  score <- Matrix::rowMeans(ga)                 # per-gene MEAN activity across cells
  out <- data.frame(gene = rownames(ga), score = as.numeric(score),
                    stringsAsFactors = FALSE)
  out <- out[out$score > 0, ]
  write.table(out, opt$`out-scores`, sep = "\t", quote = FALSE, row.names = FALSE)
  list(status = "success", n_genes = nrow(out), n_cells = length(barcodes))
}, error = function(e) list(status = "error", reason = conditionMessage(e)))

if (!is.null(opt$`out-json`)) write_json(res, opt$`out-json`)
cat(sprintf('{"status":"%s"%s}\n', res$status,
            if (!is.null(res$n_genes)) sprintf(',"n_genes":%d', res$n_genes) else
            if (!is.null(res$reason)) sprintf(',"reason":"%s"', gsub('"', "'", res$reason)) else ""))
