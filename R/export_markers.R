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
