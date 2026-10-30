log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(DESeq2); library(data.table); library(ggplot2); library(BiocParallel) })
register(MulticoreParam(snakemake@threads))

I <- snakemake@input; O <- snakemake@output
meta <- fread(I$metadata, colClasses = "character")
d <- fread(I$counts); counts <- as.matrix(d[, -1]); rownames(counts) <- d[[1]]; counts <- counts[, meta$sample]
cls <- fread(I$classification, colClasses = "character")
lnc_genes <- unique(cls$gene_id)

cd <- data.frame(row.names = meta$sample, accession = factor(meta$accession), condition = relevel(factor(meta$condition), "control"),
                 tissue = factor(meta$tissue), study = factor(meta$study))
keep <- rowSums(counts >= 10) >= 3
counts <- counts[keep, ]
covars <- c(if (nlevels(cd$tissue) > 1) "tissue", if (nlevels(cd$study) > 1 && !all(table(cd$study, cd$accession) > 0)) NULL)
base <- paste(c(covars, "accession", "condition"), collapse = " + ")

# Condition main effect
dds <- DESeqDataSetFromMatrix(counts, cd, as.formula(paste("~", base)))
dds <- DESeq(dds, parallel = TRUE)
res_c <- as.data.table(results(dds, name = resultsNames(dds)[grep("^condition_", resultsNames(dds))[1]]), keep.rownames = "feature")
res_c[, is_lncrna := feature %in% lnc_genes]
fwrite(res_c[order(pvalue)], O$condition, sep = "\t")

# G x E: LRT of full (with accession:condition) vs reduced
full <- as.formula(paste("~", base, "+ accession:condition"))
dds2 <- DESeqDataSetFromMatrix(counts, cd, full)
dds2 <- DESeq(dds2, test = "LRT", reduced = as.formula(paste("~", base)), parallel = TRUE)
res_g <- as.data.table(results(dds2), keep.rownames = "feature")
res_g[, is_lncrna := feature %in% lnc_genes]
fwrite(res_g[order(pvalue)], O$gxe, sep = "\t")

pdf(O$pdf, width = 8, height = 6)
res_c[, class := ifelse(is_lncrna, "lncRNA", "coding")]
print(ggplot(res_c[!is.na(padj)], aes(log2FoldChange, -log10(pvalue), colour = padj < 0.05)) + geom_point(alpha = .4, size = .7) +
        facet_wrap(~class) + theme_bw() + labs(title = "Condition effect"))
res_g[, class := ifelse(is_lncrna, "lncRNA", "coding")]
print(ggplot(res_g[!is.na(padj)], aes(class, fill = padj < 0.05)) + geom_bar(position = "fill") + theme_bw() +
        labs(title = "Fraction of features with significant genotype x condition interaction", y = "fraction"))
dev.off()
cat(sprintf("condition DE (padj<0.05): %d (%d lncRNA); GxE: %d (%d lncRNA)\n",
            sum(res_c$padj < 0.05, na.rm = TRUE), sum(res_c$padj < 0.05 & res_c$is_lncrna, na.rm = TRUE),
            sum(res_g$padj < 0.05, na.rm = TRUE), sum(res_g$padj < 0.05 & res_g$is_lncrna, na.rm = TRUE)))
