library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

cds <- as.cell_data_set(Merge269)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

Merge269[["Monocle_clusters"]] <- cds@clusters@listData[["UMAP"]][["clusters"]]


marker_test_res <- top_markers(cds, group_cells_by = "cluster", 
                               reference_cells=1000, cores=8)


top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(2, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

rowData(cds)$gene_short_name <- row.names(rowData(cds))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    norm_method = "log",
                    ordering_type="maximal_on_diag",
                    axis_order = "group_marker",
                    max.size=4)
