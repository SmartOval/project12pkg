library(scRNAseq)
library(SingleCellExperiment)
library(usethis)
library(ggplot2)

sce <- sce <- BaronPancreasData("mouse")
example_sce <- sce
counts_mat <- counts(sce)
cell_type <- colData(sce)$label
use_data(example_sce, overwrite = TRUE)
