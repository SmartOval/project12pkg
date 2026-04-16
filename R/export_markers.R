#' Exports all markers
#'
#' @param all_markers A DataFrame of all cells including genes, and cell-type.
#' @param marker_summary A dataframe with all distinct cell types with their total markers and top marker.
#' @param output_dir Creates a directory for the data to be exported in
#' @return An exported table of all the genes and their markers.
#' @export
export_markers <- function(all_markers, marker_summary, output_dir) {
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  utils::write.table(all_markers, file.path(output_dir, "marker_genes.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(marker_summary, file.path(output_dir, "marker_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  invisible(output_dir)
}
