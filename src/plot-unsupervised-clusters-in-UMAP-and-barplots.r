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

# group.colors <- metacell.palette[match(sapply(strsplit(metacell.labels, "_"), function(x){ paste(x[1], x[2], sep = "_") }), names(metacell.palette))]

obj_integr$seurat_clusters_annotated_sorted <- obj_integr$seurat_clusters_annotated

my_levels <- c("KC premit_c33", "KC postmit_c13", "KC postmit_c14", "KC postmit_c36",
                "Melanocyte_c0",  "Melanocyte_c1", "Melanocyte_c4",  "Melanocyte_c6", "Melanocyte_c7",
                "Melanocyte_c16", "Melanocyte_c18", "Melanocyte_c20", "Melanocyte_c23", "Melanocyte_c26", "Melanocyte_c37", "Melanocyte_c38",
                "Fb 2_c9", "Fb 2_c10", "Fb 2_c17", "Pericyte_c19", "VE 3_c12", "VE 3_c21", "LE 1_c29",
                "Tc_c5", "Th_c2", "Th_c11", "Th_c22", "Th_c30", "Th_c32", "Th_c35", "Mast cell_c3",
                "Mast cell_c25", "Mast cell_c27", "Mast cell_c31", "Plasma cell_c24", "Mac 2_c15", "DC2_c8",
                "DC2_c28", "pDC_c34" )
obj_integr$seurat_clusters_annotated_sorted <- factor(obj_integr$seurat_clusters_annotated_sorted, levels = my_levels )

# save(metacell.data, metacell.labels, metacell.colors, metacell.palette, file = "output/metacell-indata-opossom.RData")

p1 <- DimPlot(obj_integr, reduction = "umap", label = FALSE, raster = FALSE, repel = TRUE, cols = metacell.palette, group.by = "seurat_clusters_annotated_sorted") + theme(legend.position = "none")
p2 <- ggplot(obj_integr@meta.data, aes(y = seurat_clusters_annotated_sorted, fill = seurat_clusters_annotated)) +
  geom_bar() +
  theme_minimal() +
  scale_fill_manual(values = metacell.palette) +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  ) +
  xlab("Cluster") +
  ylab("Frequency")

p3 <- ggplot(obj_integr@meta.data, aes(x = factor(1), fill = seurat_clusters_annotated_sorted)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  coord_polar(theta = "y") + 
  facet_wrap(orig.ident ~ .) +
  scale_fill_manual(values = metacell.palette) +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  ) +
  xlab("Celltype") +
  ylab("Frequency")


p <- (p1 + p2) / p3 + plot_layout(nrow=2,heights=c(2,1))

DimPlot(obj_integr, reduction = "umap", label = TRUE, raster = FALSE, repel = TRUE, cols = metacell.palette, group.by = "seurat_clusters_annotated", split.by = ) + theme(legend.position = "none")
DimPlot(obj_integr, reduction = "umap", label = FALSE, raster = FALSE, repel = TRUE, cols = metacell.palette, group.by = "seurat_clusters_annotated_sorted", split.by = "orig.ident", ncol = 4) + theme(legend.position = "none")
png("output/unsupervised-clusters-in-UMAP-and-barplots.png", units = "in", width = 12, height = 8, res = 400)
p + plot_annotation(
  title = 'Graph-based clustering',
  tag_levels = 'A',
  caption = 'A: UMAP colored by seurat clusters, B: Distribution of cells among seurat clusters, C: Distribution of cells among seurat clusters per experiment'
)
dev.off()
