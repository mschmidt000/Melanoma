library(here)
library(heatmap.plus)
library(Seurat)
source(here("src","paths.r"))
source(here("src","opossom-functions.r"))

load("output/metacell-indata-opossom.RData")
filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)

marker_genes_df <- read.csv("literature/Danaher-2017_tumor-infiltrating-leukocytes/40425_2017_215_MOESM1_ESM.csv", sep = ";")
marker_genes_list <- lapply(unique(marker_genes_df$Cell.Type), function(x){ marker_genes_df$Gene[marker_genes_df$Cell.Type %in% x]})
names(marker_genes_list) <- unique(marker_genes_df$Cell.Type)

filename <- here(output_data_path, "integrated-seurat-obj.RData")
save(obj_integr, file = filename)

marker_genes_df_intersect <- marker_genes_df[which(marker_genes_df$Gene %in% rownames(metacell.data)),]
data_hm <- metacell.data[marker_genes_df_intersect$Gene,]
data_hm_gs <- t(do.call(cbind, by(data_hm, marker_genes_df_intersect$Cell.Type, colMeans))[,unique(marker_genes_df$Cell.Type)])

obj_integr <- AddModuleScore(
    obj_integr,
    marker_genes_list,
    name = "danaher-",
    search = TRUE
)

colnames(obj_integr@meta.data)[grep("danaher", colnames(obj_integr@meta.data))] <- paste("danaher", names(marker_genes_list), sep = "_")

FeaturePlot(
    obj_integr,
    grep("danaher", colnames(obj_integr@meta.data), value = TRUE)[13:14],
    order = TRUE
)

heatmap.plus(
    data_hm_gs[, names(metacell.colors)],
    ColSideColors = cbind(metacell.colors, metacell.colors)
)

