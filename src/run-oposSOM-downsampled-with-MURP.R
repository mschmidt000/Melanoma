### run opossom with metacells
### 12.10.22
library(Seurat)
library(oposSOM)
library(MURP)
library(here)
source(here("src","paths.r"))
source(here("src","seurat-functions.r"))
source(here("src","opossom-functions.r"))
source(here("src","convert-huge-sparse-matrix-to-matrix.r"))

filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)

annotated_clusters <- annotate_clusters(obj_integr, "predicted.id.Cell_type_nicknames")
obj_integr$seurat_clusters_annotated <- annotated_clusters[obj_integr$seurat_clusters]
obj_integr$seurat_clusters_annotated <- paste0(obj_integr$seurat_clusters," (",obj_integr$seurat_clusters_annotated , ")")
obj_integr$seurat_clusters_annotated <- factor(obj_integr$seurat_clusters_annotated, levels = unique(obj_integr$seurat_clusters_annotated)[order(unique(obj_integr$seurat_clusters))])

my_data <- as_matrix(GetAssayData(obj_integr, assay = "RNA", slot = "data"))
result = MURP(Data =my_data, cores = 10, iter = 1, omega = 1/20, seed = 723, fast = TRUE, cluster_iter_max = 3)
KBicPlot(murpResult = result)
MURPNestedGridPlot(murpResult = result)
