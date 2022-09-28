### run opossom on melanocyte cells
### 26.05.22
library(scrat)
library(here)
source(here("src","paths.r"))
source(here("src","seurat-functions.r"))
source(here("src","scrat-functions.r"))
source(here("src","convert-huge-sparse-matrix-to-matrix.r"))


filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)

annotated_clusters <- annotate_clusters(obj_integr, "predicted.id.Cell_type_nicknames")
obj_integr$seurat_clusters_annotated <- annotated_clusters[obj_integr$seurat_clusters]
obj_integr$seurat_clusters_annotated <- paste0(obj_integr$seurat_clusters," (",obj_integr$seurat_clusters_annotated , ")")
obj_integr$seurat_clusters_annotated <- factor(obj_integr$seurat_clusters_annotated, levels = unique(obj_integr$seurat_clusters_annotated)[order(unique(obj_integr$seurat_clusters))])

Idents(obj_integr) <- "seurat_clusters_annotated"
# obj_mel_ds <- subset(obj_mel, downsample = 300) %>%
#   ScaleData(verbose = FALSE, features = rownames(obj_mel))
# rm(obj_mel)
env <- scrat.new(list(
  dataset.name = paste(Sys.Date(), "melanoma","metacells", sep = "_"),
  dim.1stLvlSom = "auto",
  dim.2ndLvlSom = "auto",
  activated.modules = list(
    "reporting" = TRUE,
    "primary.analysis" = TRUE,
    "sample.similarity.analysis" = FALSE,
    "geneset.analysis" = TRUE,
    "psf.analysis" = TRUE,
    "group.analysis" = TRUE,
    "difference.analysis" = TRUE
  ),
  preprocessing = list(
    cellcycle.correction = FALSE,
    create.meta.cells = TRUE,
    feature.centralization = TRUE,
    sample.quantile.normalization = TRUE
  ),
  standard.spot.modules = "group.overexpression"
))


# definition of indata, group.labels and group.colors
DefaultAssay(obj_integr) <- "RNA"
var_genes_20k <- VariableFeatures(FindVariableFeatures(obj_integr, selection.method = "vst", nfeatures = 20000))
env$indata <- GetAssayData(obj_integr, slot = "counts", assay = "RNA")[var_genes_20k,] %>%
  as_matrix()
env$group.labels <- Idents(obj_integr)
env$group.colors <- env$group.labels %>%
  dplyr::n_distinct() %>%
  scales::hue_pal()(.)
names(env$group.colors) <- levels(env$group.labels) 
env$group.colors <- env$group.colors[match(env$group.labels, names(env$group.colors))]
names(env$group.colors) <- names(env$group.labels)
scrat.run(env)