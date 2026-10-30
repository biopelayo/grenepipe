log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(tximport); library(data.table) })

I <- snakemake@input; O <- snakemake@output
samples <- fread(I$samples, colClasses = "character")
tx2gene <- fread(I$tx2gene, header = FALSE, col.names = c("TXNAME", "GENEID"))
files <- setNames(unlist(I$quants), samples$sample)
stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
counts <- round(txi$counts); tpm <- txi$abundance

cls <- fread(I$classification, colClasses = "character")
lnc_genes <- unique(cls$gene_id)
coding_genes <- readLines(I$coding)
# Reference lncRNA-biotype genes are neither in coding list nor novel; keep them in "all" only.
is_lnc <- rownames(counts) %in% lnc_genes
is_cod <- rownames(counts) %in% coding_genes

write_mat <- function(m, path) fwrite(data.table(feature = rownames(m), m), path, sep = "\t")
write_mat(counts[is_lnc, , drop = FALSE], O$counts_lnc)
write_mat(tpm[is_lnc, , drop = FALSE], O$tpm_lnc)
write_mat(counts[is_cod, , drop = FALSE], O$counts_cod)
write_mat(tpm[is_cod, , drop = FALSE], O$tpm_cod)
write_mat(counts, O$counts_all)
write_mat(tpm, O$tpm_all)

meta <- samples[, .(sample, accession, accession_id, tissue, condition, study)]
meta[, is_lncrna_feature_count := sum(is_lnc)]
fwrite(meta, O$metadata, sep = "\t")
cat(sprintf("genes: %d total, %d lncRNA, %d coding; samples: %d\n", nrow(counts), sum(is_lnc), sum(is_cod), ncol(counts)))
