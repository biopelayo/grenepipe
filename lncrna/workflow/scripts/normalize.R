log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(DESeq2); library(limma); library(data.table); library(ggplot2) })

I <- snakemake@input; O <- snakemake@output; P <- snakemake@params
meta <- fread(I$metadata, colClasses = "character")
read_mat <- function(p) { d <- fread(p); m <- as.matrix(d[, -1]); rownames(m) <- d[[1]]; m[, meta$sample, drop = FALSE] }

process <- function(counts, tpm, label) {
  expressed <- rowMeans(tpm >= P$expressed_tpm) >= P$expressed_frac
  counts <- counts[expressed, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(counts, data.frame(row.names = meta$sample, study = factor(meta$study)), ~1)
  vst_m <- assay(vst(dds, blind = TRUE))
  # Batch (study) and tissue removed for unsupervised analyses only; models downstream use raw counts + covariates.
  batch <- if (length(unique(meta$study)) > 1) factor(meta$study) else NULL
  tissue <- if (length(unique(meta$tissue)) > 1) factor(meta$tissue) else NULL
  bc <- vst_m
  if (!is.null(batch) || !is.null(tissue)) {
    bc <- removeBatchEffect(vst_m, batch = batch, batch2 = tissue)
  }
  list(vst = vst_m, bc = bc, stats = data.table(set = label, features_total = nrow(tpm), features_expressed = sum(expressed)))
}
lnc <- process(read_mat(I$counts_lnc), read_mat(I$tpm_lnc), "lncRNA")
cod <- process(read_mat(I$counts_cod), read_mat(I$tpm_cod), "coding")

wm <- function(m, p) fwrite(data.table(feature = rownames(m), m), p, sep = "\t")
wm(lnc$vst, O$vst_lnc); wm(cod$vst, O$vst_cod); wm(lnc$bc, O$vst_lnc_bc); wm(cod$bc, O$vst_cod_bc)
fwrite(rbind(lnc$stats, cod$stats), O$stats, sep = "\t")

pdf(O$pca, width = 8, height = 6)
for (nm in c("lncRNA", "coding")) {
  m <- if (nm == "lncRNA") lnc$bc else cod$bc
  v <- order(matrixStats::rowVars(m), decreasing = TRUE)[1:min(2000, nrow(m))]
  pc <- prcomp(t(m[v, ]), scale. = TRUE)
  d <- data.table(PC1 = pc$x[, 1], PC2 = pc$x[, 2], study = meta$study, tissue = meta$tissue, condition = meta$condition)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  print(ggplot(d, aes(PC1, PC2, colour = study, shape = tissue)) + geom_point(size = 2) + theme_bw() +
          labs(title = paste(nm, "(batch-corrected VST, top variable features)"),
               x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")))
}
dev.off()
