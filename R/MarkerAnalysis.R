# --- Filter genes detected in < 5% of cells ---

filter_genes <- function(sce, min_detect_rate = 0.05){
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }
  if(!is.numeric(min_detect_rate) || length(min_detect_rate) != 1 ||
     min_detect_rate < 0 || min_detect_rate > 1) {
    stop("`min_detect_rate` must be a single number between 0 and 1")
  }
  counts_mat <- SingleCellExperiment::counts(sce)
  detect_rate <- rowSums(counts_mat > 0) / ncol(counts_mat)
  keep_genes <- detect_rate >= min_detect_rate
  sce[keep_genes, ]
}




detect_rate_all <- rowSums(counts_mat > 0) / ncol(counts_mat)
keep_genes <- detect_rate_all >= 0.05
counts_mat <- counts_mat[keep_genes, ]
cat("Genes after filtering:", nrow(counts_mat), "\n")




# --- Log-normalize ---



lib_sizes <- colSums(counts_mat)
mat_norm <- log2(t(t(counts_mat) / lib_sizes * 1e6) + 1)





# --- One-vs-rest Wilcoxon testing per cell type ---



cell_types <- unique(cell_type)
n_markers <- 5

marker_list <- lapply(cell_types, function(ct) {
  is_target <- cell_type == ct
  # Log2 fold change: mean(target) - mean(rest)
  mean_target <- rowMeans(mat_norm[, is_target, drop = FALSE])
  mean_rest   <- rowMeans(mat_norm[, !is_target, drop = FALSE])
  log2fc <- mean_target - mean_rest

  # Detection rate in target vs rest
  detect_target <- rowSums(counts_mat[, is_target, drop = FALSE] > 0) / sum(is_target)
  detect_rest   <- rowSums(counts_mat[, !is_target, drop = FALSE] > 0) / sum(!is_target)

  # Wilcoxon rank-sum test (subset of genes for speed)
  pvals <- sapply(seq_len(nrow(mat_norm)), function(i) {
    wilcox.test(mat_norm[i, is_target], mat_norm[i, !is_target],
                alternative = "greater")$p.value
  })

  data.frame(
    gene         = rownames(mat_norm),
    cell_type    = ct,
    log2fc       = log2fc,
    detect_target = detect_target,
    detect_rest  = detect_rest,
    pvalue       = pvals,
    padj         = p.adjust(pvals, method = "BH"),
    stringsAsFactors = FALSE
  )
})
all_markers <- do.call(rbind, marker_list)



# --- Select top N markers per cell type ---




top_markers <- do.call(rbind, lapply(cell_types, function(ct) {
  sub <- all_markers[all_markers$cell_type == ct, ]
  sub <- sub[order(sub$log2fc, decreasing = TRUE), ]
  head(sub, n_markers)
}))
