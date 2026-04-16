test_that("normalize_counts returns matrix with correct dimensions", {
  data("example_sce", package = "project12pkg")

  sce_filt <- filter_genes(example_sce)
  mat_norm <- normalize_counts(sce_filt)

  expect_true(is.matrix(mat_norm) || inherits(mat_norm, "Matrix"))
  expect_equal(dim(mat_norm), dim(SingleCellExperiment::counts(sce_filt)))
})

test_that("find markers errors on invalid group column", {
  data("example_sce", package = "project12pkg")

  sce_filt <- filter_genes(example_sce)
  mat_norm <- normalize_counts(sce_filt)

  expect_error(
    find_markers(sce_filt, mat_norm, group_col = "not_a_column"),
    "`group_col` must be a valid column in colData\\(sce\\)."
  )
})

