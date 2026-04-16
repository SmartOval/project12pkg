test_that("filter_genes returns a SingleCellExperiment Object", {
  data("example_sce", package = "project12pkg")

  out <- filter_genes(example_sce, min_detect_rate = 0.05)

  expect_true(inherits(out, "SingleCellExperiment"))
})

test_that("filter_genes errors on invalid sce input", {
  expect_error(
    filter_genes(5),
    "`sce` must be a SingleCellExperiment object."
  )
})

test_that("marker pipline returns a non-empty marker table", {
  data("example_sce", package = "project12pkg")

  sce_filt <- filter_genes(example_sce, min_detect_rate = 0.05)
  mat_norm <- normalize_counts(sce_filt)
  all_markers <- find_markers(sce_filt, mat_norm, group_col = "label")

  expect_s3_class(all_markers, "data.frame")
  expect_true(nrow(all_markers) > 0)
  expect_true(all(c("gene", "cell_type", "log2fc", "pvalue", "padj") %in% colnames(all_markers)))
})
