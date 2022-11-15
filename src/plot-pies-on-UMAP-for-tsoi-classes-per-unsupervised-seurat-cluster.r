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

pdf("output/pies-for-tsoi-class-per-unsupervised-melanocyte-cluster.pdf")
ggplot(mels_sub@meta.data, aes(x = seurat_clusters_annotated_sorted, fill = Tsoi_classification)) +
  geom_bar(position = "fill") + 
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_brewer(palette="Set1") +
  facet_wrap( ~seurat_clusters_annotated)
dev.off()

