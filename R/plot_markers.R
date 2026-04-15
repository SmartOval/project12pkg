build_dotplot <- function(all_markers, top_markers, mat_norm, cell_type) {
  if (!is.data.frame(all_markers) || !is.data.frame(top_markers)) {
    stop("`all_markers` and `top_markers` must be data.frames.")
  }
  gene_levels <- unique(as.character(top_markers$gene))
  ct_levels   <- as.character(cell_types)
  plot_df <- all_markers[all_markers$gene %in% gene_levels, ]

  # Add mean expression per (gene, cell_type) combination
  plot_df$mean_expr <- mapply(function(g, ct) {
    mean(as.numeric(mat_norm[g, cell_type == ct]))
  }, plot_df$gene, plot_df$cell_type)

  plot_df$gene      <- factor(plot_df$gene,      levels = rev(gene_levels))
  plot_df$cell_type <- factor(plot_df$cell_type, levels = ct_levels)
  plot_df
}



plot_markers <- function(plot_df) {
  if (!is.data.frame(plot_df)) {
    stop("`plot_df` must be a data.Frame.")
  }
  ggplot2:ggplot(plot_df, aes(x = cell_type, y = gene,
                      size = detect_target, color = mean_expr)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(range = c(1, 6), name = "Detection Rate") +
    ggplot2::scale_color_viridis_c(name = "Mean Expression") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "Top Marker Genes per Cell Type")
}
