### reduce dimension and find clusters
### 19.05.22

library(clustree)

n_dims_use <- 30

filename <- here(figures_path, "02_dimension-reduction-and-clustering.pdf")
pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)

for(i in 1:length(runs_list)){
  perplexity <- sqrt(ncol(runs_list[[runs[i]]]@assays$RNA@counts))
  runs_list[[runs[i]]] <- runs_list[[runs[i]]] %>%
                                RunPCA( features = VariableFeatures(runs_list[[runs[i]]])) %>%
                                FindNeighbors(dims = 1:n_dims_use) %>%
                                FindNeighbors(dims = 1:n_dims_use) %>%
                                FindClusters(resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1)) %>%
                                RunTSNE(dims = 1:n_dims_use, perplexity = perplexity) %>%
                                RunUMAP(dims = 1:n_dims_use )  %>%
                                CellCycleScoring(
                                  s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes,
                                  set.ident = FALSE
                                )

  plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))
  mtext(runs[i], side=3, line = -2, cex=3, at=-0.04, font=3,adj=0)
  
  plot(VizDimLoadings(runs_list[[runs[i]]], dims = 1:2, reduction = "pca"))
  plot(DimPlot(runs_list[[runs[i]]], reduction = "pca"))
  plot(ElbowPlot(runs_list[[runs[i]]]))
  plot(FeaturePlot(runs_list[[runs[i]]], reduction = "tsne", 
                   features = c("percent.mt", "percent.rps", "nCount_RNA", "nFeature_RNA" )))
  plot(DimPlot(runs_list[[runs[i]]], reduction = "tsne", label = TRUE, group.by = c("seurat_clusters", "Phase")))
  plot(DimPlot(runs_list[[runs[i]]], reduction = "umap", label = TRUE, group.by = c("seurat_clusters", "Phase")))
  clustree( x = runs_list[[runs[i]]]@meta.data, prefix = "RNA_snn_res.") + ggtitle("Clustering Tree visualizing Clusterings at Different Clustering-Resolutions")
 }

dev.off()


for(i in seq_along(runs_list)){
  seurat_obj <- runs_list[[runs[i]]]
  filename <- here(output_data_path,paste0(runs[i], "_seurat-obj.RData"))
  save(seurat_obj, file = filename)
}
rm(seurat_obj)
