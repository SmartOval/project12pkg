#!/usr/bin/env Rapp
#| name: project12pkg
#| title: Marker gene CLI for Project 12
#| description: Filter genes, identify marker genes, summarize markers, and export results.

suppressPackageStartupMessages({
  library(project12pkg)
  library(utils)
  library(stats)
  library(SummarizedExperiment)
})

switch(
  "",

  #| title: Run marker gene analysis
  #| description: Run the marker pipeline on the bundled example dataset and export results.
  run = {
    #| description: Minimum detection rate for gene filtering
    #| short: d
    min_detect_rate <- 0.05

    #| description: Number of top markers per cell type
    #| short: n
    n_markers <- 5L

    #| description: Output directory
    #| short: o
    output_dir <- ""

    if (output_dir == "") {
      stop("--output-dir is required", call. = FALSE)
    }

    if (min_detect_rate < 0 || min_detect_rate > 1) {
      stop("--min-detect-rate must be between 0 and 1", call. = FALSE)
    }

    if (n_markers < 1) {
      stop("--n-markers must be at least 1", call. = FALSE)
    }

    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    data("example_sce", package = "project12pkg", envir = environment())

    sce_filt <- filter_genes(example_sce, min_detect_rate = min_detect_rate)
    mat_norm <- normalize_counts(sce_filt)
    all_markers <- find_markers(sce_filt, mat_norm, group_col = "label")
    top_markers <- select_top_markers(all_markers, n_markers = n_markers)
    marker_summary <- summarize_markers(all_markers)

    export_markers(all_markers, marker_summary, output_dir)

    utils::write.table(
      top_markers,
      file.path(output_dir, "top_markers.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    message("Done.")
    message("Wrote:")
    message("  ", file.path(output_dir, "marker_genes.tsv"))
    message("  ", file.path(output_dir, "marker_summary.tsv"))
    message("  ", file.path(output_dir, "top_markers.tsv"))
  },

  #| title: Validate inputs
  #| description: Check that CLI arguments are valid before running analysis.
  validate = {
    #| description: Minimum detection rate for gene filtering
    #| short: d
    min_detect_rate <- 0.05

    #| description: Number of top markers per cell type
    #| short: n
    n_markers <- 5L

    #| description: Output directory
    #| short: o
    output_dir <- ""

    if (output_dir == "") {
      stop("--output-dir is required", call. = FALSE)
    }

    if (min_detect_rate < 0 || min_detect_rate > 1) {
      stop("--min-detect-rate must be between 0 and 1", call. = FALSE)
    }

    if (n_markers < 1) {
      stop("--n-markers must be at least 1", call. = FALSE)
    }

    message("Inputs look valid.")
    message("min_detect_rate = ", min_detect_rate)
    message("n_markers = ", n_markers)
    message("output_dir = ", output_dir)
  }
)
