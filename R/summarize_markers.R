summarize_markers <- function(all_markers, padj_cutoff = 1, log2fc_cutoff = 1) {
  if (!is.data.frame(all_markers)) {
    stop("`all_markers` must be a data.Frame.")
  }
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
