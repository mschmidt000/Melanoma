### subclustering all melanocytes
### 20.05.22

Idents(obj_integr) <- "predicted.id.Cell_type_nicknames"
obj_mel <- subset(obj_integr, idents = "Melanocyte")
Idents(obj_integr) <- "orig.ident"
DefaultAssay(obj_mel) <- "RNA"
integr_features <- rownames(obj_mel)
runs_list_mel <- SplitObject(obj_mel, split.by = "orig.ident")
reference <- which(names(runs_list_mel) %in% c("LE-583-KM", "LE-585-HM"))
integr_anchors <- FindIntegrationAnchors(runs_list_mel, dims = 1:30, reference = reference, anchor.features = 10000)
rm(runs_list_mel)
obj_mel <- IntegrateData(anchorset = integr_anchors, dims = 1:30)
DefaultAssay(obj_mel) <- "integrated"

obj_mel <- FindVariableFeatures(obj_mel)
obj_mel <- ScaleData(obj_mel, verbose = FALSE, features = rownames(obj_mel)) %>%
  RunPCA(npcs = n_dims_use, verbose = FALSE) %>%
  RunTSNE(reduction = "pca", dims = 1:n_dims_use, perplexity = sqrt(ncol(obj_mel))) %>%
  RunUMAP(dims = 1:n_dims_use) %>%
  FindNeighbors(reduction = "pca", dims = 1:n_dims_use) %>%
  FindClusters(resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2))

filename <- here(output_data_path, "integrated-seurat-obj-melanocytes.RData")
save(obj_mel, file = filename)

filename <- here(figures_path, "05_dimension-reduction-and-clustering-of-mels.pdf")
pdf(filename, height = 7, width = 14)
DimPlot(obj_mel, reduction = "umap", label = TRUE, group.by = "orig.ident", repel = TRUE)
barplot(table(obj_mel@meta.data$orig.ident), main = paste0("Number of Cells per Sample"), col = scales::hue_pal()(length(table(obj_mel@meta.data$orig.ident))))
DimPlot(obj_mel, reduction = "umap", label = TRUE, group.by = paste0("integrated_snn_res.", c(0.5, 0.6, 0.8, 1)))
clustree::clustree(x = obj_mel@meta.data, prefix = "integrated_snn_res.")
DimPlot(obj_mel, reduction = "umap", label = FALSE, repel = TRUE, group.by = "Phase")
FeaturePlot(obj_mel, reduction = "umap", order = TRUE, features = c("percent.mt", "percent.rps", "nCount_RNA", "nFeature_RNA"))
dev.off()

interesting_idents_clust_mels <- paste("predicted.id.Cell_type_nicknames", "subcl_mels", sep = "_")

obj_integr@meta.data[, interesting_idents_clust_mels] <- as.character(obj_integr@meta.data[, "predicted.id.Cell_type_nicknames"])
# obj_integr@meta.data[colnames(obj_mel), interesting_idents_clust_mels] <- paste(as.character(obj_mel@meta.data[, interesting_idents_long[1]]), obj_mel@meta.data[, "integrated_snn_res.1"], sep = "_")
melanocytes <- grep("Melanocyte", obj_integr@meta.data[, "predicted.id.Cell_type_nicknames"])
obj_integr@meta.data[melanocytes, interesting_idents_clust_mels] <- paste("Mel", obj_integr@meta.data[melanocytes, "integrated_snn_res.1.2"], sep = "_")
obj_integr@meta.data[, interesting_idents_clust_mels] <- factor(obj_integr@meta.data[, interesting_idents_clust_mels], 
                                                                levels = c(levels(obj_integr$predicted.id.Cell_type_nicknames)[1:4], unique(obj_integr@meta.data[melanocytes, interesting_idents_clust_mels]),levels(obj_integr$predicted.id.Cell_type_nicknames)[6:39])
                                                                  )
filename <- here(output_data_path, "integrated-seurat-obj.RData")
save(obj_integr, file = filename)
rm(obj_mel)

Idents(obj_integr) <- interesting_idents_clust_mels

filename <- here(figures_path, "05_integrated-dimension-reduction-and-clustering-with-mel-subclusters.pdf")
pdf(filename, height = 7, width = 14)
DimPlot(obj_integr, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()
