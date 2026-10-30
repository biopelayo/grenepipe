log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(MatrixEQTL); library(data.table); library(rtracklayer) })

I <- snakemake@input; O <- snakemake@output; E <- snakemake@params$eqtl
meta <- fread(I$metadata, colClasses = "character")
d <- fread(I$vst); X <- as.matrix(d[, -1]); rownames(X) <- d[[1]]; X <- X[, meta$sample]

# accession-level expression (mean over samples of the same accession_id)
ok <- meta$accession_id != ""
Xa <- sapply(split(which(ok), meta$accession_id[ok]), function(idx) rowMeans(X[, idx, drop = FALSE]))

gt <- fread(I$gt); G <- as.matrix(gt[, -1]); rownames(G) <- gt[[1]]
common <- intersect(colnames(Xa), colnames(G))
cat(sprintf("accessions with expression and genotypes: %d\n", length(common)))
if (length(common) < 20) stop("fewer than 20 accessions with both expression and genotypes")
Xa <- Xa[, common]; G <- G[, common]
G <- G[rowMeans(is.na(G)) <= 0.1, ]
G[is.na(G)] <- rowMeans(G, na.rm = TRUE)[row(G)[is.na(G)]]

# hidden factors: top PCs of expression (PEER-like)
pcs <- prcomp(t(Xa), scale. = TRUE)$x[, seq_len(min(E$n_hidden_factors, length(common) - 2)), drop = FALSE]
cvrt <- SlicedData$new(); cvrt$CreateFromMatrix(t(pcs))

# gene positions: TSS +/- window from lncRNA GTF (gene-level)
gtf <- as.data.table(import(I$gtf)); gtf <- gtf[type %in% c("transcript", "exon")]
genepos <- gtf[, .(chr = as.character(seqnames)[1], s = min(start), e = max(end)), by = gene_id]
genepos <- genepos[gene_id %in% rownames(Xa)]
Xa <- Xa[genepos$gene_id, ]
pos <- fread(I$pos); pos <- pos[snp %in% rownames(G)]
G <- G[pos$snp, ]

expr <- SlicedData$new(); expr$CreateFromMatrix(Xa); expr$ResliceCombined(2000)
snps <- SlicedData$new(); snps$CreateFromMatrix(G); snps$ResliceCombined(5000)

me <- Matrix_eQTL_main(snps = snps, gene = expr, cvrt = cvrt, output_file_name = NULL, output_file_name.cis = NULL,
                       pvOutputThreshold = 0, pvOutputThreshold.cis = 1e-3, useModel = modelLINEAR,
                       snpspos = as.data.frame(pos[, .(snp, chr, pos)]),
                       genepos = as.data.frame(genepos[, .(gene_id, chr, s, e)]),
                       cisDist = E$cis_window, verbose = FALSE, pvalue.hist = FALSE, min.pv.by.genesnp = TRUE)
cis <- as.data.table(me$cis$eqtls)
setnames(cis, c("snps", "gene"), c("snp", "gene_id"))
fwrite(cis[order(pvalue)], O$cis, sep = "\t")

# eGenes: per-gene min p, Bonferroni within gene by number of cis SNPs tested, then BH across genes
ntest <- as.data.table(me$cis$min.pv.gene, keep.rownames = "gene_id")
setnames(ntest, 2, "min_p")
nsnp <- cis[, .N, by = gene_id]  # approximate count (SNPs passing output threshold); conservative fallback below
genes_tested <- data.table(gene_id = names(me$cis$min.pv.gene), min_p = as.numeric(me$cis$min.pv.gene))
n_cis <- sapply(genes_tested$gene_id, function(g) { r <- genepos[gene_id == g]; sum(pos$chr == r$chr & pos$pos >= r$s - E$cis_window & pos$pos <= r$e + E$cis_window) })
genes_tested[, p_bonf := pmin(1, min_p * pmax(1, n_cis))]
genes_tested[, fdr := p.adjust(p_bonf, "BH")]
egenes <- genes_tested[fdr < E$fdr][order(fdr)]
egenes <- merge(egenes, cis[, .SD[which.min(pvalue)], by = gene_id][, .(gene_id, lead_snp = snp, beta, lead_p = pvalue)], by = "gene_id")
fwrite(egenes, O$egenes, sep = "\t")

pdf(O$pdf, width = 8, height = 6)
hist(cis$pvalue, breaks = 50, main = "cis p-value distribution (p < 1e-3 reported)", xlab = "p")
if (nrow(egenes)) {
  plot(density(abs(egenes$beta)), main = sprintf("|beta| of lead cis-eQTL, %d lncRNA eGenes (FDR < %.2f)", nrow(egenes), E$fdr))
}
dev.off()
cat(sprintf("cis tests: %d genes, %d SNPs; eGenes: %d\n", nrow(Xa), nrow(G), nrow(egenes)))
