log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(ConsensusClusterPlus); library(data.table); library(mclust); library(cluster) })

I <- snakemake@input; O <- snakemake@output; P <- snakemake@params; C <- P$consensus
meta <- fread(I$metadata, colClasses = "character")
read_mat <- function(p) { d <- fread(p); m <- as.matrix(d[, -1]); rownames(m) <- d[[1]]; m[, meta$sample, drop = FALSE] }

# Stratification is at the accession level: average replicates / tissues per accession.
collapse <- function(m) {
  acc <- meta$accession
  s <- sapply(split(seq_along(acc), acc), function(idx) rowMeans(m[, idx, drop = FALSE]))
  s
}
top_var <- function(m, n) { v <- matrixStats::rowVars(m); m[order(v, decreasing = TRUE)[1:min(n, nrow(m))], ] }

# PAC (proportion of ambiguous clustering) to pick k, Senbabaoglu et al. 2014
pac <- function(cc, k, lo = 0.1, hi = 0.9) {
  M <- cc[[k]]$consensusMatrix; x <- M[lower.tri(M)]
  mean(x > lo & x < hi)
}

run_cc <- function(m, label) {
  m <- sweep(m, 1, rowMeans(m))
  cc <- ConsensusClusterPlus(m, maxK = C$max_k, reps = C$reps, pItem = C$p_item, pFeature = 1,
                             clusterAlg = C$cluster_alg, distance = C$distance, seed = C$seed, plot = NULL)
  summ <- rbindlist(lapply(2:C$max_k, function(k) {
    cl <- cc[[k]]$consensusClass
    sil <- mean(silhouette(cl, as.dist(1 - cor(m)))[, 3])
    data.table(set = label, k = k, PAC = pac(cc, k), mean_silhouette = sil,
               cluster_sizes = paste(table(cl), collapse = ";"))
  }))
  list(cc = cc, summary = summ)
}

lnc <- top_var(collapse(read_mat(I$vst_lnc)), P$top_variable)
cod <- top_var(collapse(read_mat(I$vst_cod)), P$top_variable)
r_lnc <- run_cc(lnc, "lncRNA")
r_cod <- run_cc(cod, "coding")
summary <- rbind(r_lnc$summary, r_cod$summary)
best_k <- function(s) s[which.min(s$PAC), k]
k_lnc <- best_k(r_lnc$summary); k_cod <- best_k(r_cod$summary)
summary[, selected := (set == "lncRNA" & k == k_lnc) | (set == "coding" & k == k_cod)]
fwrite(summary, O$summary, sep = "\t")

cl_lnc <- r_lnc$cc[[k_lnc]]$consensusClass
cl_cod <- r_cod$cc[[k_cod]]$consensusClass
clusters <- data.table(accession = names(cl_lnc), lncrna_cluster = cl_lnc, coding_cluster = cl_cod[names(cl_lnc)])
acc_meta <- unique(meta[, .(accession, accession_id)])
clusters <- merge(clusters, acc_meta, by = "accession", all.x = TRUE)

cmp <- data.table(comparison = "lncRNA_vs_coding", k_lncrna = k_lnc, k_coding = k_cod,
                  ARI = adjustedRandIndex(cl_lnc, cl_cod[names(cl_lnc)]))
if (length(I$labels)) {
  lab <- fread(I$labels, colClasses = "character")
  clusters <- merge(clusters, lab, by = "accession_id", all.x = TRUE)
  ok <- !is.na(clusters$label)
  cmp <- rbind(cmp, data.table(comparison = "lncRNA_vs_external_label", k_lncrna = k_lnc, k_coding = NA,
                               ARI = adjustedRandIndex(clusters$lncrna_cluster[ok], clusters$label[ok])),
               data.table(comparison = "coding_vs_external_label", k_lncrna = NA, k_coding = k_cod,
                          ARI = adjustedRandIndex(clusters$coding_cluster[ok], clusters$label[ok])))
}
fwrite(clusters, O$clusters, sep = "\t")
fwrite(cmp, O$compare, sep = "\t")

pdf(O$pdf, width = 9, height = 7)
for (r in list(r_lnc, r_cod)) {
  plot(2:C$max_k, r$summary$PAC, type = "b", xlab = "k", ylab = "PAC", main = paste(r$summary$set[1], "- PAC by k"))
  for (k in 2:C$max_k) {
    M <- r$cc[[k]]$consensusMatrix; o <- order(r$cc[[k]]$consensusClass)
    image(M[o, o], main = paste(r$summary$set[1], "consensus matrix k =", k), axes = FALSE, col = hcl.colors(50, "Blues", rev = TRUE))
  }
}
dev.off()
cat(sprintf("lncRNA: k=%d (PAC %.3f); coding: k=%d; ARI=%.3f\n", k_lnc, min(r_lnc$summary$PAC), k_cod, cmp$ARI[1]))
