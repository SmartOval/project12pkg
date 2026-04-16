library(scRNAseq)
library(SingleCellExperiment)
library(usethis)
library(ggplot2)
sce <- sce <- BaronPancreasData("mouse")
example_sce <- sce

data("example_sce", package = "project12pkg")

sce_filt <- filter_genes(example_sce)
mat_norm <- normalize_counts(sce_filt)
all_markers <- find_markers(sce_filt, mat_norm)

top_markers <- select_top_markers(all_markers, n_markers = 5)
marker_summary <- summarize_markers(all_markers)

cell_type <- SummarizedExperiment::colData(sce_filt)[["label"]]
plot_df <- build_dotplot(all_markers, top_markers, mat_norm, cell_type)
p <- plot_markers(plot_df)

export_markers(all_markers, marker_summary, tempdir())
