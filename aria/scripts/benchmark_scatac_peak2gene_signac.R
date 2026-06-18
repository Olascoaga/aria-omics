#!/usr/bin/env Rscript
# scATAC P4.2 — Signac LinkPeaks reference driver for the peak-to-gene concordance
# benchmark. Runs in aria-bench-signac-env. Consumes the IDENTICAL paired matrices
# the export stage wrote (atac_counts.mtx / rna_counts.mtx + peaks/genes/barcodes),
# builds a Seurat object with an RNA assay and a Signac ChromatinAssay, sets the gene
# annotation from the local GTF, computes per-peak GC via the UCSC hg38 BSgenome,
# runs LinkPeaks, and writes the reference links TSV (gene, peak, score, zscore,
# pvalue) that the pure-numpy P4.2 scorer consumes.
#
# No fabrication: any missing dependency / asset is an honest non-success status in
# the output JSON, never an invented link set (mirrors the ADR-002 discipline of the
# SnapATAC2 driver).

suppressWarnings(suppressMessages({
  ok <- requireNamespace("optparse", quietly = TRUE)
}))
if (!ok) { cat('{"status":"not_run","reason":"optparse unavailable"}\n'); quit(status = 0) }
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--in-dir", type = "character"),
  make_option("--gtf", type = "character"),
  make_option("--out-links", type = "character"),
  make_option("--out-json", type = "character"),
  make_option("--distance", type = "double", default = 5e5),
  make_option("--pvalue-cutoff", type = "double", default = 0.05),
  make_option("--score-cutoff", type = "double", default = 0.0),
  make_option("--genome", type = "character", default = "BSgenome.Hsapiens.UCSC.hg38")
)))

write_json <- function(obj, path) {
  # tiny dependency-free JSON writer (flat object of scalars/strings)
  parts <- vapply(names(obj), function(k) {
    v <- obj[[k]]
    if (is.numeric(v)) sprintf('"%s": %s', k, format(v, scientific = FALSE))
    else sprintf('"%s": "%s"', k, gsub('"', "'", as.character(v)))
  }, character(1))
  writeLines(paste0("{\n  ", paste(parts, collapse = ",\n  "), "\n}"), path)
}

fail <- function(reason) {
  if (!is.null(opt$`out-json`)) write_json(list(status = "not_run", reason = reason), opt$`out-json`)
  cat(sprintf('{"status":"not_run","reason":"%s"}\n', reason))
  quit(status = 0)
}

need <- c("Matrix", "Seurat", "Signac", "GenomicRanges", "rtracklayer", "IRanges", "S4Vectors")
for (pkg in need) if (!requireNamespace(pkg, quietly = TRUE)) fail(paste0(pkg, " unavailable"))
if (!requireNamespace(opt$genome, quietly = TRUE)) fail(paste0(opt$genome, " unavailable"))

suppressWarnings(suppressMessages({
  library(Matrix); library(Seurat); library(Signac)
  library(GenomicRanges); library(rtracklayer)
}))

indir <- opt$`in-dir`
res <- tryCatch({
  barcodes <- readLines(file.path(indir, "barcodes.tsv"))
  genes <- readLines(file.path(indir, "genes.tsv"))
  peaks_df <- read.table(file.path(indir, "peaks.tsv"), sep = "\t",
                         col.names = c("chrom", "start", "end", "name"),
                         stringsAsFactors = FALSE)
  peak_ids <- paste0(peaks_df$chrom, "-", peaks_df$start, "-", peaks_df$end)

  atac <- readMM(file.path(indir, "atac_counts.mtx"))  # peaks x cells
  rna <- readMM(file.path(indir, "rna_counts.mtx"))    # genes x cells
  rownames(atac) <- peak_ids; colnames(atac) <- barcodes
  rownames(rna) <- genes;     colnames(rna) <- barcodes

  obj <- CreateSeuratObject(counts = as(rna, "CsparseMatrix"), assay = "RNA")
  obj <- NormalizeData(obj, verbose = FALSE)

  ranges <- GRanges(seqnames = peaks_df$chrom,
                    ranges = IRanges(start = peaks_df$start, end = peaks_df$end))
  chrom_assay <- CreateChromatinAssay(counts = as(atac, "CsparseMatrix"),
                                      ranges = ranges, sep = c("-", "-"))
  obj[["peaks"]] <- chrom_assay
  DefaultAssay(obj) <- "peaks"
  obj <- RunTFIDF(obj, verbose = FALSE)

  # Gene annotation from the local GTF (matches the peak chr* naming).
  gtf <- rtracklayer::import(opt$gtf)
  annot <- gtf[gtf$type == "gene"]
  if (!"gene_biotype" %in% colnames(mcols(annot)) && "gene_type" %in% colnames(mcols(annot)))
    annot$gene_biotype <- annot$gene_type
  if (!"tx_id" %in% colnames(mcols(annot))) annot$tx_id <- annot$gene_id
  annot$type <- "exon"
  Annotation(obj[["peaks"]]) <- annot

  genome <- getExportedValue(opt$genome, opt$genome)
  obj <- RegionStats(obj, genome = genome, assay = "peaks")

  use_genes <- intersect(unique(as.character(annot$gene_name)), genes)
  obj <- LinkPeaks(obj, peak.assay = "peaks", expression.assay = "RNA",
                   genes.use = use_genes, distance = opt$distance,
                   pvalue_cutoff = opt$`pvalue-cutoff`, score_cutoff = opt$`score-cutoff`,
                   verbose = FALSE)

  links <- Links(obj[["peaks"]])
  ld <- as.data.frame(links)
  out <- data.frame(gene = ld$gene, peak = ld$peak, score = ld$score,
                    zscore = ld$zscore, pvalue = ld$pvalue, stringsAsFactors = FALSE)
  write.table(out, opt$`out-links`, sep = "\t", quote = FALSE, row.names = FALSE)
  list(status = "success", n_links = nrow(out),
       n_genes = length(unique(out$gene)), n_cells = length(barcodes),
       n_peaks = length(peak_ids))
}, error = function(e) list(status = "error", reason = conditionMessage(e)))

if (!is.null(opt$`out-json`)) write_json(res, opt$`out-json`)
cat(sprintf('{"status":"%s"%s}\n', res$status,
            if (!is.null(res$n_links)) sprintf(',"n_links":%d', res$n_links) else
            if (!is.null(res$reason)) sprintf(',"reason":"%s"', gsub('"', "'", res$reason)) else ""))
