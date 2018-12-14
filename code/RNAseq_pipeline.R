## select top g genes from a gene expression matrix X
## rows of X: genes
## columns of X: samples
top_genes_index = function (g, X) {
  return(order(rowSums(X), decreasing = TRUE)[1 : g])
}

## convert read counts to lcpm from a RNA-seq read count matrix r
## rows of r: genes
## columns of r: samples
lcpm = function (r) {
  R = colSums(r)
  t(log2(((t(r) + 0.5) / (R + 1)) * 10^6))
}

## Function to implement a widely-used RNA-seq count data analysis pipeline
## for gene differential expression studies in a A/B testing setting
## Input a count matrix and a group allocation
## this pipeline gives, for each gene j
## an estimate of log2-fold differential expression X_j
## a p-value testing for differential expression
## a z-score which is N(0, 1) distributed if no differential expression
## an effective standard deviation s_j
count_to_summary = function (counts, group) {
  design = model.matrix(~group)
  dgecounts = edgeR::calcNormFactors(edgeR::DGEList(counts = counts, group = design[, 2]))
  v = limma::voom(dgecounts, design, plot = FALSE)
  lim = limma::lmFit(v)
  r.ebayes = limma::eBayes(lim)
  p = r.ebayes$p.value[, 2]
  t = r.ebayes$t[, 2]
  z = -sign(t) * qnorm(p / 2)
  X = lim$coefficients[, 2]
  s = X / z
  return (list(X = X, s = s, zscore = z, pval = p))
}
