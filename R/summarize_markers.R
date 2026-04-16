#' Summary: overlap counts between groups
#'
#' @param all_markers A DataFrame of all cells including genes, and cell-type.
#' @param padj_cutoff A value for the desired cut off for the p-adjusted value
#' @param log2fc_cutoff A value for the desired cut off for the log2fc value.
#' @return A dataframe with all distinct cell types with their total markers and top marker.
#' @export
summarize_markers <- function(all_markers, padj_cutoff = 1, log2fc_cutoff = 1) {
  if (!is.data.frame(all_markers)) {
    stop("`all_markers` must be a data.Frame.")
  }
  cell_types <- unique(all_markers$cell_type)
  data.frame(
    cell_type   = cell_types,
    n_sig       = sapply(cell_types, function(ct) {
      sum(all_markers$cell_type == ct & all_markers$padj < padj_cutoff &
            all_markers$log2fc > log2fc_cutoff)
    }),
    top_marker  = sapply(cell_types, function(ct) {
      sub <- all_markers[all_markers$cell_type == ct, ]
      sub$gene[which.max(sub$log2fc)]
    }),
    stringsAsFactors = FALSE
  )
}
