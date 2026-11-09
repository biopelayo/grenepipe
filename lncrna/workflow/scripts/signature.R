log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(glmnet); library(data.table); library(pheatmap) })

I <- snakemake@input; O <- snakemake@output; S <- snakemake@params$sig
set.seed(snakemake@params$seed)

d <- fread(I$vst_lnc); X <- as.matrix(d[, -1]); rownames(X) <- d[[1]]
clusters <- fread(I$clusters, colClasses = c(accession = "character"))
# collapse samples to accessions as in stratify.R (mean per accession)
md <- fread(I$metadata, colClasses = "character")
X <- X[, md$sample, drop = FALSE]
Xa <- sapply(split(seq_len(ncol(X)), md$accession), function(idx) rowMeans(X[, idx, drop = FALSE]))
Xa <- t(Xa[, clusters$accession, drop = FALSE])  # accessions x features
y <- factor(clusters$lncrna_cluster)
if (nlevels(y) < 2 || min(table(y)) < 3) stop("need >= 2 clusters with >= 3 accessions each for a signature")

fam <- if (nlevels(y) == 2) "binomial" else "multinomial"
nonzero_features <- function(cv, s = "lambda.1se") {
  cf <- coef(cv, s = s)
  if (!is.list(cf)) cf <- list(cf)
  feats <- unique(unlist(lapply(cf, function(b) { b <- as.matrix(b); rownames(b)[-1][b[-1, 1] != 0] })))
  feats
}
outer <- sample(rep(seq_len(S$outer_folds), length.out = nrow(Xa)))
perf <- rbindlist(lapply(seq_len(S$outer_folds), function(f) {
  tr <- outer != f; te <- !tr
  cv <- cv.glmnet(Xa[tr, ], y[tr], family = fam, alpha = S$alpha, nfolds = S$inner_folds, type.measure = "class")
  pred <- predict(cv, Xa[te, , drop = FALSE], s = "lambda.1se", type = "class")
  data.table(fold = f, n_test = sum(te), accuracy = mean(pred == as.character(y[te])),
             n_features = length(nonzero_features(cv)))
}))
perf <- rbind(perf, data.table(fold = NA, n_test = nrow(Xa), accuracy = mean(perf$accuracy), n_features = mean(perf$n_features)))
fwrite(perf, O$performance, sep = "\t")

final <- cv.glmnet(Xa, y, family = fam, alpha = S$alpha, nfolds = S$inner_folds, type.measure = "class")
cf <- coef(final, s = "lambda.1se")
if (is.list(cf)) {
  sig <- rbindlist(lapply(names(cf), function(k) { b <- as.matrix(cf[[k]]); data.table(feature = rownames(b)[-1], class = k, coef = b[-1, 1]) }))
} else {
  b <- as.matrix(cf); sig <- data.table(feature = rownames(b)[-1], class = levels(y)[2], coef = b[-1, 1])
}
sig <- sig[coef != 0]
sig <- sig[, .(max_abs_coef = max(abs(coef)), classes = paste(class, collapse = ";"), coefs = paste(round(coef, 3), collapse = ";")), by = feature][order(-max_abs_coef)]
cls <- fread(I$classification, colClasses = "character")
sig <- merge(sig, unique(cls[, .(feature = gene_id, lncrna_type, te_like = if ("te_like" %in% names(cls)) te_like else NA)]), by = "feature", all.x = TRUE)[order(-max_abs_coef)]
fwrite(sig, O$signature, sep = "\t")

pdf(O$pdf, width = 10, height = 8)
if (nrow(sig) >= 2) {
  m <- t(scale(Xa[, sig$feature, drop = FALSE]))
  ann <- data.frame(cluster = y, row.names = rownames(Xa))
  pheatmap(m, annotation_col = ann, show_colnames = nrow(Xa) <= 60, main = sprintf("lncRNA signature (%d features), accessions", nrow(sig)))
}
dev.off()
cat(sprintf("signature: %d lncRNAs, nested-CV accuracy %.3f\n", nrow(sig), perf[is.na(fold), accuracy]))
