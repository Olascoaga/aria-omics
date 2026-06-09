#!/usr/bin/env Rscript
# One-time bootstrap: materialize a SEQC/MAQC A1 reference bundle from the
# `seqc` Bioconductor data package (counts + TaqMan qPCR truth).
#
# Runs inside aria-bench-env. Writes counts.tsv / samples.tsv / taqman.tsv to
# the output dir; the Python wrapper adds manifest.json with provenance hashes.
#
# Source: seqc package (SEQC/MAQC-III, Nature Biotechnology 2014; TaqMan from
# MAQC-I, Nature Biotechnology 2006). Counts: one site's RefSeq gene table
# (default BGI), lane/flow-cell columns summed to the 4 library replicates per
# sample A/B/C/D, ERCC spike-in rows dropped. TaqMan: log2(mean A / mean B) of
# the per-replicate POG values, by gene Symbol. ARIA fabricates nothing here.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "."
count_table <- if (length(args) >= 2) args[[2]] else "ILM_refseq_gene_BGI"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(seqc))

# --- counts ---------------------------------------------------------------
data(list = count_table, package = "seqc")
x <- get(count_table)

# Split off the ERCC spike-in rows (written separately for the dose-response
# lane); the gene matrix excludes them.
ercc <- if (!is.null(x$IsERCC)) x[x$IsERCC, , drop = FALSE] else x[0, ]
if (!is.null(x$IsERCC)) x <- x[!x$IsERCC, , drop = FALSE]

cnt_cols <- grep("^[ABCD]_[0-9]+_L", colnames(x), value = TRUE)
rep_id <- sub("^([ABCD]_[0-9]+)_.*$", "\\1", cnt_cols)
reps <- unique(rep_id)
mat <- vapply(
  reps,
  function(r) rowSums(x[, cnt_cols[rep_id == r], drop = FALSE]),
  numeric(nrow(x))
)
storage.mode(mat) <- "integer"
counts <- data.frame(gene = as.character(x$Symbol), mat, check.names = FALSE)
write.table(counts, file.path(out_dir, "counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- samples --------------------------------------------------------------
samples <- data.frame(sample = reps, group = sub("_.*$", "", reps))
write.table(samples, file.path(out_dir, "samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- ERCC spike-in counts (same sample-replicate aggregation) -------------
if (nrow(ercc) > 0) {
  emat <- vapply(
    reps,
    function(r) rowSums(ercc[, cnt_cols[rep_id == r], drop = FALSE]),
    numeric(nrow(ercc))
  )
  storage.mode(emat) <- "integer"
  ercc_counts <- data.frame(ercc_id = as.character(ercc$Symbol), emat,
                            check.names = FALSE)
  write.table(ercc_counts, file.path(out_dir, "ercc_counts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# --- taqman truth (log2 A/B) ----------------------------------------------
data(taqman, package = "seqc")
a_cols <- grep("^A[0-9]+_value$", colnames(taqman), value = TRUE)
b_cols <- grep("^B[0-9]+_value$", colnames(taqman), value = TRUE)
a_mean <- rowMeans(taqman[, a_cols, drop = FALSE], na.rm = TRUE)
b_mean <- rowMeans(taqman[, b_cols, drop = FALSE], na.rm = TRUE)
log2_ab <- log2(a_mean / b_mean)
tq <- data.frame(gene = as.character(taqman$Symbol), log2_ab = log2_ab)
tq <- tq[is.finite(tq$log2_ab), , drop = FALSE]
write.table(tq, file.path(out_dir, "taqman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf(
  "SEQC_BUNDLE_DONE genes=%d samples=%d taqman=%d ercc=%d table=%s\n",
  nrow(counts), nrow(samples), nrow(tq), nrow(ercc), count_table
))
