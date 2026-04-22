#' Builds a data frame of all genes and their markers to be plotted

#' @param all_markers A DataFrame of all cells including genes, and cell-type.
#' @param top_markers A DataFrame of the top markers per gene
#' @param mat_norm Large dgeMatrix of the normalized counts of the SingleCellExperiment object
#' @param cell_type A dataframe of the filtered genes and all of their respective cell types
#' @return A dataframe of each filtered gene, ordered by cell type.
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_se)
#'
#' # We should assume preprocessing functions were run earlier
#' # all_markers <- find_all_markers(example_se)
#' # top_markers <- get_top_markers(all_markers)
#' # mat_norm <- get_normalized_matrix(example_se)
#' # cell_type <- colData(example_se)$cell_type
#'
#' # Build dotplot dataframe
#' plot_df <- build_dotplot(
#'   all_markers = all_markers,
#'   top_markers = top_markers,
#'   mat_norm = mat_norm,
#'   cell_type = cell_type
#' )
#' }
#' @export
build_dotplot <- function(all_markers, top_markers, mat_norm, cell_type) {
  if (!is.data.frame(all_markers) || !is.data.frame(top_markers)) {
    stop("`all_markers` and `top_markers` must be data.frames.")
  }
  gene_levels <- unique(as.character(top_markers$gene))
  ct_levels   <- unique(as.character(cell_type))
  plot_df <- all_markers[all_markers$gene %in% gene_levels, ]

  # Add mean expression per (gene, cell_type) combination
  plot_df$mean_expr <- mapply(function(g, ct) {
    mean(as.numeric(mat_norm[g, cell_type == ct]))
  }, plot_df$gene, plot_df$cell_type)

  plot_df$gene      <- factor(plot_df$gene,      levels = rev(gene_levels))
  plot_df$cell_type <- factor(plot_df$cell_type, levels = ct_levels)
  plot_df
}


#' Plots genes and their markers
#'
#' @param plot_df A dataframe of all filtered genes
#' @return A ggplot dotplot of the plot_df data frame
#' @examples
#' \dontrun{
#' data(example_se)
#'
#' # Example workflow
#' # plot_df <- build_dotplot(...)
#'
#' # Plot markers
#' plot_markers(plot_df)
#' }
#' @importFrom rlang .data
#' @export
plot_markers <- function(plot_df) {
  if (!is.data.frame(plot_df)) {
    stop("`plot_df` must be a data.frame.")
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["cell_type"]],
      y = .data[["gene"]],
      size = .data[["detect_target"]],
      color = .data[["mean_expr"]]
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(range = c(1, 6), name = "Detection Rate") +
    ggplot2::scale_color_viridis_c(name = "Mean Expression") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(title = "Top Marker Genes per Cell Type")
}
