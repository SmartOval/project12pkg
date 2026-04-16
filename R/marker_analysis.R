#' Filter genes by detection rate
#'
#' @param sce A SingleCellExperiment object.
#' @param min_detect_rate Minimum fraction of cells in which a gene must be detected.
#'
#' @return A filtered SingleCellExperiment object.
#' @export
filter_genes <- function(sce, min_detect_rate = 0.05){
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }
  if(!is.numeric(min_detect_rate) || length(min_detect_rate) != 1 ||
     min_detect_rate < 0 || min_detect_rate > 1) {
    stop("`min_detect_rate` must be a single number between 0 and 1")
  }
  counts_mat <- SingleCellExperiment::counts(sce)
  detect_rate <- Matrix::rowSums(counts_mat > 0) / ncol(counts_mat)
  keep_genes <- detect_rate >= min_detect_rate
  sce[keep_genes, ]
}

#' Normalize gene counts
#'
#' @param sce A SingleCellExperiment object.
#'
#' @return A normalized SingleCellExperiment object.
#' @export
normalize_counts <- function(sce) {
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }

  counts_mat <- SingleCellExperiment::counts(sce)
  lib_sizes <- Matrix::colSums(counts_mat)

  if (any(lib_sizes == 0)) {
    stop("One or more cells have library size 0.")
  }
  log2(Matrix::t(Matrix::t(counts_mat) / lib_sizes * 1e6) + 1)
}

#' Filter genes by detection rate
#'
#' @param sce A SingleCellExperiment object.
#' @param mat_norm A normalized expression matrix.
#' @param group_col Column name in \code{colData(sce)} containing cell type labels.
#' @return A filtered SingleCellExperiment object.
#' @export
find_markers <- function(sce, mat_norm, group_col = "label") {
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }
  if (!(is.matrix(mat_norm) || inherits(mat_norm, "Matrix"))) {
    stop("`mat_norm` must be a matrix-like object.")
  }
  if(!group_col %in% colnames(SummarizedExperiment::colData(sce))) {
    stop("`group_col` must be a valid column in colData(sce).")
  }
  cell_type <- SummarizedExperiment::colData(sce)[[group_col]]
  counts_mat <- SingleCellExperiment::counts(sce)
  cell_types <- unique(cell_type)
  marker_list <- lapply(cell_types, function(ct) {
    is_target <- cell_type == ct
    # Log2 fold change: mean(target) - mean(rest)
    mean_target <- Matrix::rowMeans(mat_norm[, is_target, drop = FALSE])
    mean_rest   <- Matrix::rowMeans(mat_norm[, !is_target, drop = FALSE])
    log2fc <- mean_target - mean_rest

    # Detection rate in target vs rest
    detect_target <- Matrix::rowSums(counts_mat[, is_target, drop = FALSE] > 0) / sum(is_target)
    detect_rest   <- Matrix::rowSums(counts_mat[, !is_target, drop = FALSE] > 0) / sum(!is_target)

    # Wilcoxon rank-sum test (subset of genes for speed)
    pvals <- sapply(seq_len(nrow(mat_norm)), function(i) {
      stats::wilcox.test(mat_norm[i, is_target], mat_norm[i, !is_target],
                  alternative = "greater")$p.value
    })

    data.frame(
      gene         = rownames(mat_norm),
      cell_type    = ct,
      log2fc       = log2fc,
      detect_target = detect_target,
      detect_rest  = detect_rest,
      pvalue       = pvals,
      padj         = stats::p.adjust(pvals, method = "BH"),
      stringsAsFactors = FALSE
    )
  })
  all_markers <- do.call(rbind, marker_list)
  all_markers

}


#' Select top N markers per cell type
#'
#' @param all_markers A DataFrame of all cells including genes, and cell-type.
#' @param n_markers Preferred number of markers.
#'
#' @return A cell_types object with inputted number of top markers per cell types.
#' @export
select_top_markers <- function(all_markers, n_markers = 5) {
  if (!is.data.frame(all_markers)) {
    stop("`all_makers` must be a data.frame.")
  }
  if (!is.numeric(n_markers) || length(n_markers) != 1 || n_markers < 1) {
    stop("`n_markers` must be a positive number.")
  }
  cell_types <- unique(all_markers$cell_type)
  do.call(rbind, lapply(cell_types, function(ct) {
    sub <- all_markers[all_markers$cell_type == ct, ]
    sub <- sub[order(sub$log2fc, decreasing = TRUE), ]
    utils::head(sub, n_markers)
  }))
}
