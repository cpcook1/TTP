library(Seurat)
DataMix <- readRDS("E:/Human - 10 clusters/yale-method/07062021-organ-T-integration/T_organ_integration_res1.rds")
Trm <- subset(DataMix, idents = "0")

library(monocle)
MyMTX <- Trm@assays[["RNA"]]@data
GeneAnno <- as.data.frame(row.names(MyMTX))
colnames(GeneAnno) <- "gene_short_name"
row.names(GeneAnno) <- GeneAnno$gene_short_name

MyClust <- Trm@meta.data$seurat_clusters
MyClust <- as.numeric(MyClust)
MyClust <- as.data.frame(MyClust)
colnames(MyClust) <- "clusters"
myPCA <- as.data.frame(Trm@reductions$pca@cell.embeddings)
rownames(MyClust) <- rownames(myPCA)

cell.cluster <- MyClust

Ha <- data.frame(do.call('rbind', strsplit(as.character(row.names(cell.cluster)),'_',fixed=TRUE)))[1]

clusts <- paste("cl.",as.character(cell.cluster$clusters),sep="")
cell.cluster <- cbind(cell.cluster,Ha,clusts)
colnames(cell.cluster) <- c("Clusts","iCellR.Conds","iCellR.Clusts")
Samp <- new("AnnotatedDataFrame", data = cell.cluster)
Anno <- new("AnnotatedDataFrame", data = GeneAnno)
my.monoc.obj <- newCellDataSet(as.matrix(MyMTX),phenoData = Samp, featureData = Anno)

## find disperesedgenes 
my.monoc.obj <- estimateSizeFactors(my.monoc.obj)
my.monoc.obj <- estimateDispersions(my.monoc.obj)
my.monoc.obj <- detectGenes(my.monoc.obj, min_expr = 3 )
disp_table <- dispersionTable(my.monoc.obj)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
my.monoc.obj <- setOrderingFilter(my.monoc.obj, unsup_clustering_genes$gene_id)

# tSNE
my.monoc.obj <- reduceDimension(my.monoc.obj, max_components = 2, num_dim = 5,reduction_method = 'tSNE', verbose = T)
# cluster 
my.monoc.obj <- clusterCells(my.monoc.obj, num_clusters = 5)

## plot conditions and clusters based on iCellR analysis 
A <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "iCellR.Conds")
B <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "iCellR.Clusts")

## plot clusters based monocle analysis 
C <- plot_cell_clusters(my.monoc.obj, 1, 2, color = "Cluster")

# get marker genes from iCellR analysis
library(dplyr)
print(head(fData(my.monoc.obj)))
expressed_genes <- row.names(subset(fData(my.monoc.obj), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(my.monoc.obj[expressed_genes,])
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
my.monoc.obj <- setOrderingFilter(my.monoc.obj, ordering_genes)
plot_ordering_genes(my.monoc.obj)

my.monoc.obj <- reduceDimension(my.monoc.obj, max_components = 2, method = 'DDRTree')
# order cells 
my.monoc.obj <- orderCells(my.monoc.obj)

# plot based on iCellR analysis and marker genes from iCellR
plot_cell_trajectory(my.monoc.obj, color_by = "iCellR.Clusts")
plot_cell_trajectory(my.monoc.obj, color_by = "iCellR.Conds")
plot_cell_trajectory(my.monoc.obj, color_by = "State")
plot_cell_trajectory(my.monoc.obj, color_by = "Pseudotime")
plot_cell_trajectory(my.monoc.obj, color_by = "Pseudotime")
## heatmap genes from iCellR

plot_pseudotime_heatmap(my.monoc.obj[ordering_genes,],
                        num_clusters = 3,
                        cores = 1,
                        cluster_rows = T,
                        use_gene_short_name = T,
                        show_rownames = T)

gene1 <- my.monoc.obj@featureData@data[["gene_short_name"]] %in% c("CD8A", "CD4", "GZMA")
plot_genes_jitter(my.monoc.obj[gene1,],
                  grouping = "State",
                  min_expr = 0.1)

plot_genes_in_pseudotime(my.monoc.obj[gene1,], color_by = "State")
