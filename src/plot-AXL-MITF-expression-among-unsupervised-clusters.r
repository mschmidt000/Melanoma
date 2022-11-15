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

obj_integr$Mel_seurat <- rep(FALSE, ncol(obj_integr))
obj_integr$Mel_seurat[grep("Melanocyte", obj_integr$seurat_clusters_annotated)] <- TRUE
obj_integr$Mel_seurat[grep("Mast", obj_integr$seurat_clusters_annotated)] <- TRUE
table(obj_integr$Mel_seurat)
mels_sub <- subset(obj_integr, Mel_seurat)

pdf("output/AXL-MITF-expression-among-unsupervised-clusters.pdf", height = 4)
DotPlot(mels_sub, features = c("AXL", "MITF"), group.by = "seurat_clusters_annotated") +
  RotatedAxis() +
  coord_flip() +
  scale_color_gradient2(high = "red", low = "blue") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
dev.off()

