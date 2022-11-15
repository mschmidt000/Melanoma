library(Seurat)
library(here)
library(ggplot2)
library(patchwork)
source(here("src","paths.r"))
source(here("src","seurat-functions.r"))
source(here("src","opossom-functions.r"))

filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)

mels_sub <- subset(obj_integr, Mel)

FeaturePlot(
    mels_sub,
    grep("Tsoi", colnames(obj_integr@meta.data), value = TRUE)[1:7],
    order = TRUE
)


