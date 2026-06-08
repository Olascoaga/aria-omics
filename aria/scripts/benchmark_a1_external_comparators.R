#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--counts", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--output-dir", type = "character", dest = "output_dir"),
  make_option("--output-json", type = "character", dest = "output_json"),
  make_option("--padj", type = "double", default = 0.05),
  make_option("--lfc", type = "double", default = 0.5)
)
opt <- parse_args(OptionParser(option_list = option_list))

write_payload <- function(payload) {
  write_json(payload, opt$output_json, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

fail <- function(error_type, details) {
  write_payload(list(status = "error", error_type = error_type, details = details))
  quit(status = 0)
}

required <- c("counts", "metadata", "output_dir", "output_json")
missing_args <- required[vapply(required, function(x) is.null(opt[[x]]) || opt[[x]] == "", logical(1))]
if (length(missing_args) > 0) {
  fail("InvalidArguments", paste("Missing required argument(s):", paste(missing_args, collapse = ", ")))
}

missing_packages <- setdiff(
  c("DESeq2", "edgeR", "limma", "apeglm"),
  rownames(installed.packages())
)
if (length(missing_packages) > 0) {
  fail("MissingDependency", paste("Missing R package(s):", paste(missing_packages, collapse = ", ")))
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

counts_dt <- read.delim(opt$counts, check.names = FALSE)
meta_dt <- read.delim(opt$metadata, check.names = FALSE)
if (!"gene" %in% names(counts_dt)) {
  fail("InvalidInput", "Counts table must contain a 'gene' column.")
}
if (!all(c("sample", "condition") %in% names(meta_dt))) {
  fail("InvalidInput", "Metadata table must contain 'sample' and 'condition' columns.")
}

genes <- counts_dt$gene
count_mat <- as.matrix(counts_dt[, setdiff(names(counts_dt), "gene"), drop = FALSE])
rownames(count_mat) <- genes
mode(count_mat) <- "integer"

meta <- as.data.frame(meta_dt)
rownames(meta) <- meta$sample
meta <- meta[colnames(count_mat), , drop = FALSE]
if (any(is.na(meta$condition))) {
  fail("InvalidInput", "Metadata samples do not match count-matrix columns.")
}
meta$condition <- relevel(factor(meta$condition), ref = "COND_A")

standardize <- function(df, method, out_path) {
  df <- as.data.frame(df)
  names(df) <- make.names(names(df))
  if (!"gene" %in% names(df)) {
    df$gene <- rownames(df)
  }
  if (!"log2FoldChange" %in% names(df) && "logFC" %in% names(df)) {
    df$log2FoldChange <- df$logFC
  }
  if (!"padj" %in% names(df) && "FDR" %in% names(df)) {
    df$padj <- df$FDR
  }
  if (!"padj" %in% names(df) && "adj.P.Val" %in% names(df)) {
    df$padj <- df$adj.P.Val
  }
  if (!"pvalue" %in% names(df) && "PValue" %in% names(df)) {
    df$pvalue <- df$PValue
  }
  if (!"pvalue" %in% names(df) && "P.Value" %in% names(df)) {
    df$pvalue <- df$P.Value
  }
  keep <- intersect(c("gene", "log2FoldChange", "pvalue", "padj"), names(df))
  df <- df[, keep, drop = FALSE]
  write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  list(
    status = "success",
    method = method,
    results_tsv = normalizePath(out_path, mustWork = FALSE),
    n_tested = nrow(df),
    n_called = sum(!is.na(df$padj) & df$padj < opt$padj &
                     !is.na(df$log2FoldChange) & abs(df$log2FoldChange) >= opt$lfc)
  )
}

methods <- list()

tryCatch({
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = meta, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "COND_B", "COND_A")))
  res$gene <- rownames(res)
  methods$deseq2 <- standardize(res, "DESeq2", file.path(opt$output_dir, "deseq2.tsv"))
}, error = function(e) {
  methods$deseq2 <<- list(status = "error", method = "DESeq2", error_type = class(e)[1], details = conditionMessage(e))
})

tryCatch({
  suppressPackageStartupMessages(library(edgeR))
  group <- meta$condition
  design <- model.matrix(~ group)
  y <- DGEList(counts = count_mat, group = group)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tab <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none")$table)
  tab$gene <- rownames(tab)
  methods$edgeR_QLF <- standardize(tab, "edgeR-QLF", file.path(opt$output_dir, "edgeR_QLF.tsv"))
}, error = function(e) {
  methods$edgeR_QLF <<- list(status = "error", method = "edgeR-QLF", error_type = class(e)[1], details = conditionMessage(e))
})

tryCatch({
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  group <- meta$condition
  design <- model.matrix(~ group)
  y <- DGEList(counts = count_mat, group = group)
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tab <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tab$gene <- rownames(tab)
  methods$limma_voom <- standardize(tab, "limma-voom", file.path(opt$output_dir, "limma_voom.tsv"))
}, error = function(e) {
  methods$limma_voom <<- list(status = "error", method = "limma-voom", error_type = class(e)[1], details = conditionMessage(e))
})

write_payload(list(
  status = "success",
  methods = methods,
  thresholds = list(padj = opt$padj, lfc = opt$lfc),
  versions = list(
    R = paste(R.version$major, R.version$minor, sep = "."),
    DESeq2 = as.character(packageVersion("DESeq2")),
    edgeR = as.character(packageVersion("edgeR")),
    limma = as.character(packageVersion("limma"))
  )
))
