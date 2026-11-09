log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(ConsensusClusterPlus); library(data.table); library(mclust); library(cluster) })

# Accession stratification on three representations (docs/04_logica_comparativa.md, block B3):
#   coding_vst   : coding genes, batch-corrected VST, Pearson      (baseline: genetic/tissue structure)
#   lncrna_vst   : lncRNAs, batch-corrected VST, Pearson
#   lncrna_bin   : lncRNAs, ON/OFF (TPM > on_tpm) per accession, Jaccard  (the biology: most lncRNAs are OFF)
# plus ARI against each other and against every column of an external accession metadata table.

I <- snakemake@input; O <- snakemake@output; P <- snakemake@params; C <- P$consensus
ON_TPM <- if (!is.null(P$on_tpm)) P$on_tpm else 0.5
meta <- fread(I$metadata, colClasses = "character")
read_mat <- function(p) { d <- fread(p); m <- as.matrix(d[, -1]); rownames(m) <- d[[1]]; m[, meta$sample, drop = FALSE] }
collapse <- function(m) sapply(split(seq_along(meta$accession), meta$accession), function(idx) rowMeans(m[, idx, drop = FALSE]))
top_var <- function(m, n) { v <- matrixStats::rowVars(m); m[order(v, decreasing = TRUE)[1:min(n, nrow(m))], , drop = FALSE] }
pac <- function(cc, k, lo = 0.1, hi = 0.9) { M <- cc[[k]]$consensusMatrix; x <- M[lower.tri(M)]; mean(x > lo & x < hi) }
jaccard_dist <- function(B) {            # B: loci x accessions, logical
  B <- B * 1; inter <- crossprod(B); tot <- outer(colSums(B), colSums(B), "+") - inter
  D <- 1 - inter / pmax(tot, 1); diag(D) <- 0; as.dist(D)
}

run_cc <- function(x, label, dist_obj = NULL) {
  if (is.null(dist_obj)) {
    x <- sweep(x, 1, rowMeans(x))
    cc <- ConsensusClusterPlus(x, maxK = C$max_k, reps = C$reps, pItem = C$p_item, pFeature = 1,
                               clusterAlg = C$cluster_alg, distance = C$distance, seed = C$seed, plot = NULL)
    dd <- as.dist(1 - cor(x))
  } else {
    cc <- ConsensusClusterPlus(dist_obj, maxK = C$max_k, reps = C$reps, pItem = C$p_item,
                               clusterAlg = "pam", seed = C$seed, plot = NULL)
    dd <- dist_obj
  }
  summ <- rbindlist(lapply(2:C$max_k, function(k) {
    cl <- cc[[k]]$consensusClass
    data.table(set = label, k = k, PAC = pac(cc, k), mean_silhouette = mean(silhouette(cl, dd)[, 3]),
               cluster_sizes = paste(table(cl), collapse = ";"))
  }))
  list(cc = cc, summary = summ, k = summ[which.min(PAC), k])
}

# ---- matrices -----------------------------------------------------------------------------------
lnc_vst <- top_var(collapse(read_mat(I$vst_lnc)), P$top_variable)
cod_vst <- top_var(collapse(read_mat(I$vst_cod)), P$top_variable)
tpm_acc <- collapse(read_mat(I$tpm_lnc))
B <- tpm_acc > ON_TPM
n_on <- rowSums(B); n_acc <- ncol(B)
B <- B[n_on >= 2 & n_on <= n_acc - 2, , drop = FALSE]      # informative ON/OFF loci
cat(sprintf("accessions: %d; lncRNA ON/OFF informative loci: %d of %d (ON in 2..n-2)\n", n_acc, nrow(B), nrow(tpm_acc)))

runs <- list(coding_vst = run_cc(cod_vst, "coding_vst"), lncrna_vst = run_cc(lnc_vst, "lncrna_vst"))
if (nrow(B) >= 20) runs$lncrna_bin <- run_cc(NULL, "lncrna_bin", dist_obj = jaccard_dist(B))

summary <- rbindlist(lapply(runs, `[[`, "summary"))
summary[, selected := FALSE]
for (nm in names(runs)) summary[set == nm & k == runs[[nm]]$k, selected := TRUE]
fwrite(summary, O$summary, sep = "\t")

accs <- colnames(lnc_vst)
clusters <- data.table(accession = accs)
for (nm in names(runs)) clusters[[paste0(nm, "_cluster")]] <- runs[[nm]]$cc[[runs[[nm]]$k]]$consensusClass[accs]
clusters <- merge(clusters, unique(meta[, .(accession, accession_id)]), by = "accession", all.x = TRUE)

# per-accession silencing-efficiency proxy: number of lncRNA loci ON (block B5)
clusters[, n_lncrna_on := colSums(tpm_acc > ON_TPM)[accession]]

# ---- comparisons --------------------------------------------------------------------------------
cmp <- list()
pairs <- combn(names(runs), 2, simplify = FALSE)
for (p in pairs) cmp[[length(cmp) + 1]] <- data.table(a = p[1], b = p[2], k_a = runs[[p[1]]]$k, k_b = runs[[p[2]]]$k,
  ARI = adjustedRandIndex(clusters[[paste0(p[1], "_cluster")]], clusters[[paste0(p[2], "_cluster")]]), test = "")
if (length(I$labels)) {
  lab <- fread(I$labels, colClasses = "character")
  clusters <- merge(clusters, lab, by = "accession_id", all.x = TRUE)
  set.seed(C$seed)
  for (col in setdiff(names(lab), c("accession_id", "name"))) {
    v <- clusters[[col]]; ok <- !is.na(v) & v != ""
    if (sum(ok) < 10) next
    num <- suppressWarnings(as.numeric(v[ok])); is_num <- all(!is.na(num)) && length(unique(num)) > 6
    for (nm in names(runs)) {
      cl <- clusters[[paste0(nm, "_cluster")]][ok]
      if (is_num) {                                  # continuous covariate: Kruskal-Wallis across clusters
        pv <- tryCatch(kruskal.test(num, factor(cl))$p.value, error = function(e) NA)
        cmp[[length(cmp) + 1]] <- data.table(a = nm, b = col, k_a = runs[[nm]]$k, k_b = NA, ARI = NA, test = sprintf("kruskal p=%.3g", pv))
      } else {                                       # categorical label: ARI + permutation p
        ari <- adjustedRandIndex(cl, v[ok]); perm <- replicate(1000, adjustedRandIndex(cl, sample(v[ok])))
        cmp[[length(cmp) + 1]] <- data.table(a = nm, b = col, k_a = runs[[nm]]$k, k_b = length(unique(v[ok])), ARI = ari,
                                             test = sprintf("perm p=%.3g", (sum(perm >= ari) + 1) / 1001))
      }
    }
  }
}
cmp <- rbindlist(cmp)
fwrite(clusters, O$clusters, sep = "\t")
fwrite(cmp, O$compare, sep = "\t")

pdf(O$pdf, width = 9, height = 7)
for (r in runs) {
  plot(2:C$max_k, r$summary$PAC, type = "b", xlab = "k", ylab = "PAC", main = paste(r$summary$set[1], "- PAC by k"))
  for (k in 2:C$max_k) {
    M <- r$cc[[k]]$consensusMatrix; o <- order(r$cc[[k]]$consensusClass)
    image(M[o, o], main = paste(r$summary$set[1], "consensus matrix k =", k), axes = FALSE, col = hcl.colors(50, "Blues", rev = TRUE))
  }
}
if (!is.null(runs$lncrna_bin)) {
  hist(clusters$n_lncrna_on, breaks = 20, main = "lncRNA loci ON per accession (silencing-efficiency proxy)", xlab = "n ON")
}
dev.off()
for (nm in names(runs)) cat(sprintf("%s: k=%d (PAC %.3f)\n", nm, runs[[nm]]$k, min(runs[[nm]]$summary$PAC)))
print(cmp)
