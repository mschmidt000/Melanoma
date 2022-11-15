library(Seurat)
library(here)
library(ggplot2)
library(patchwork)
source(here("src","paths.r"))
source(here("src","seurat-functions.r"))
source(here("src","opossom-functions.r"))

filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)
load("output/metacell-indata-opossom.RData")

filename <- here(output_data_path, "reynolds-skin-marker.RData")
load(filename)

png("output/marker-dotplot-for-unsupervised-clusters.png", units = "in", width = 15, height = 10, res = 400)
DotPlot(obj_integr, features = skin_marker, group.by = "seurat_clusters_annotated_sorted", dot.scale = 6) +
  scale_color_gradient2(high = "red", low = "blue") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)
  )
dev.off()

