pacman::p_load(Seurat, tidyverse, here, ggplot2, patchwork)
source(here("src","seurat-functions.r"))
source(here("src","liana-functions.r"))
source(here("src","paths.r"))

filename <- file.path(output_data_path, "integrated-seurat-obj.RData")
load(filename)

obj_integr$Mel_seurat <- rep(FALSE, ncol(obj_integr))
obj_integr$Mel_seurat[grep("Melanocyte", obj_integr$seurat_clusters_annotated)] <- TRUE
obj_integr$Mel_seurat[grep("Mast", obj_integr$seurat_clusters_annotated)] <- TRUE
table(obj_integr$Mel_seurat)
mels_sub <- subset(obj_integr, Mel_seurat)

pdf("output/pies-for-cc-phases-per-unsupervised-melanocyte-cluster.pdf")
ggplot(mels_sub@meta.data, aes(x = seurat_clusters_annotated_sorted, fill = Phase)) +
  geom_bar(position = "fill") + 
  coor_polar(theta = "y")
  theme_void() +
  scale_fill_brewer(palette="Set1") +
  theme(
    # legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  ) +
  xlab("Cluster") +
  ylab("Frequency")
dev.off()

