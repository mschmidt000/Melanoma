### integrating good quality data sets
### 20.05.22

runs_to_integrate <- c( "s1_LE479_NM_T", "LE489KG_Rep", "LE493BR_Rep", "LE497NA_Rep", "LE497ST_Rep", "LE-501-DM_Rep", 
                        "LE-511-MW", "LE_569_HH", "LE_577_WR", "LE_579_DE", "LE-583-KM", "LE-585-HM" )

runs_list_subset <- sapply(runs_to_integrate, function(x){ runs_list[[x]] })
integration_features <- sapply(runs_list_subset, function(x){ rownames(x)})
integration_features <- Reduce(intersect, integration_features)
integr_anchors <- FindIntegrationAnchors(runs_list_subset, dims = 1:30, reference = c(11, 12), anchor.features = integration_features)
obj_integr <- IntegrateData(anchorset = integr_anchors, dims = 1:30)
DefaultAssay(obj_integr) <- "integrated"
filename <- file.path(output_data_path, "integrated_seuratObj.RData")
save(obj_integr, file = filename)

obj_integr <- ScaleData(obj_integr, verbose = FALSE, features = rownames(obj_integr)) %>%
                RunPCA(npcs = n_dims_use, verbose = FALSE) %>%
                RunTSNE(reduction = "pca", dims = 1:n_dims_use, perplexity = sqrt(ncol(obj_integr))) %>%
                RunUMAP(dims = 1:n_dims_use) %>%
                FindNeighbors(reduction = "pca", dims = 1:n_dims_use) %>%
                FindClusters(resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1))

filename <- file.path(output_data_path, "integrated_seuratObj.RData")
save(obj_integr, file = filename)

colors_list <- sapply( colnames(obj_integr)[grep("predicted.id.", colnames(obj_integr))], function(x){
  colors_list[["predicted.id.Cell_type"]] <- scales::hue_pal()(length(levels(obj_integr$seurat_clusters))) 
                                                                  })

colors_list[["seurat_clusters"]] <- scales::hue_pal()(length(levels(obj_integr@meta.data$seurat_clusters)))
names(colors_list[["seurat_clusters"]]) <- levels(obj_integr@meta.data$seurat_clusters)
colors_list <- c(colors_list, 

colors_list[["predicted.id.Cell_type"]] <- scales::hue_pal()(length(levels(obj_integr@meta.data$seurat_clusters)))
names(colors_list[["predicted.id.Cell_type"]]) <- levels(obj_integr@meta.data$predicted.id.Cell_type)



seurat_clusters_colors <- scales::hue_pal()(length(levels(obj_integr@meta.data$seurat_clusters)))
names(seurat_clusters_colors) <- levels(obj_integr@meta.data$seurat_clusters)
seurat_clusters_colors <- pals::polychrome(length(levels(obj_integr@meta.data$seurat_clusters)))
names(seurat_clusters_colors) <- levels(obj_integr@meta.data$seurat_clusters)


pdf("ANALYSIS_Apr22_seurat_obj.integr.11.data.sets_Analysis.pdf", width = 14)


# DimPlot(obj_integr, reduction = "pca", label = TRUE, group.by = "orig.ident", repel = TRUE) + ggtitle("PCA colored by Run")
# DimPlot(obj_integr, reduction = "tsne", label = TRUE, group.by = "orig.ident", repel = TRUE) + ggtitle("tSNE colored by Run")
DimPlot(obj_integr, reduction = "umap", label = TRUE, group.by = "orig.ident", repel = TRUE) + ggtitle("UMAP colored by Run")
DimPlot(obj_integr, reduction = "umap", label = FALSE, group.by = "orig.ident", repel = TRUE) + ggtitle("UMAP colored by Run")
barplot(table(obj_integr@meta.data$orig.ident), main = paste0("Number of Cells per Sample"), col = scales::hue_pal()(length(table(obj_integr@meta.data$orig.ident))))
# labs( caption = "PCA colored by Louvain-Clustering results based on different resolutions")
# my_plot <- DimPlot(obj_integr, reduction = "pca", label = TRUE, group.by = paste0( "integrated_snn_res.", c(0.5, 0.6, 0.8, 1) ))
# my_plot + patchwork:::plot_annotation("PCA colored by Louvain-Clustering results based on different Clustering-Resolutions (0.5, 0.6, 0.8, 1)")
# my_plot <- DimPlot(obj_integr, reduction = "tsne", label = TRUE, group.by = paste0( "integrated_snn_res.", c(0.5, 0.6, 0.8, 1) ))
# my_plot + patchwork:::plot_annotation("tSNE colored by Louvain-Clustering results based on different Clustering-Resolutions (0.5, 0.6, 0.8, 1)")
my_plot <- DimPlot(obj_integr, reduction = "umap", label = TRUE, group.by = paste0( "integrated_snn_res.", c(0.5, 0.6, 0.8, 1) ))
my_plot + patchwork:::plot_annotation("UMAP colored by Louvain-Clustering results based on different Clustering-Resolutions (0.5, 0.6, 0.8, 1)")
clustree( x = obj_integr@meta.data, prefix = "integrated_snn_res.") + ggtitle("Clustering Tree visualizing Clusterings at Different Clustering-Resolutions")
# DimPlot(obj_integr, reduction = "pca", label = TRUE, group.by = "Phase") + ggtitle("PCA colored by Cell Cycle Phase")
# DimPlot(obj_integr, reduction = "tsne", label = TRUE, group.by = "Phase") + ggtitle("tSNE colored by Cell Cycle Phase")
DimPlot(obj_integr, reduction = "umap", label = FALSE, repel = TRUE, group.by = "Phase") + ggtitle("UMAP colored by Cell Cycle Phase")
my_plot <- FeaturePlot(obj_integr, reduction = "umap", order = TRUE, features = c("percent.mt", "percent.rps", "nCount_RNA", "nFeature_RNA" ))
my_plot + patchwork:::plot_annotation("UMAP colored by the Percentage of all counts belonging to mitochondrial genes (percent.mt) or ribosomal genes (percent.rps) and colored by number of counts (nCount_RNA) or features (nFeature_RNA)",
                                      theme = theme(plot.title = element_text(size = 10)))
DotPlot(obj_integr, features = c("percent.mt", "percent.rps" ), col.min = 0, col.max = 100)
# runs <- unique(obj_integr$orig.ident)
# layout(matrix(c(1:4),ncol=2, byrow = TRUE),widths=c(1,1))
# for(run in runs){
# idys <- grep(run, obj_integr@meta.data$orig.ident )
# barplot(table(obj_integr@meta.data$seurat_cluster[idys]), main = paste0("Number of Cells per Cluster (",run,")"), col = seurat_clusters_colors)
# }
dev.off()

interesting_idents_long <- c("predicted.id.Cell_type", "predicted.id.Cell_group", "predicted.id.Flow_gate" )
interesting_idents_short <- c("Reynolds 21 Cell Type", "Reynolds 21 Cell Group", "Reynolds 21 Flow Gate" )

my_palettes <- sapply( interesting_idents_long, function(x){
  temp <- alpha(pals::polychrome(length(unique(obj_integr@meta.data[, x]))), 0.8)
  names(temp) <- sort(unique(obj_integr@meta.data[, x]))
  temp						
})

Idents(obj_integr) <- "orig.ident"
filename <- "ANALYSIS_Apr22_obj.integr.11.data.sets_predicted.ids.from.reynolds.pdf"
pdf(filename, height = 7, width = 14)
plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))
mtext("All 11 Data Sets", side=3, line = -2, cex=3, at=-0.04, font=3,adj=0)

for(i in seq_along(interesting_idents_long)){
  obj_integr@meta.data[, interesting_idents_long[i]] <- factor( obj_integr@meta.data[, interesting_idents_long[i]], level = sort(unique(obj_integr@meta.data[, interesting_idents_long[i]])))
  my_title <- paste0("UMAP colored by ", interesting_idents_short[i] )
  print(DimPlot(obj_integr, reduction = "umap", label = TRUE, group.by = interesting_idents_long[i], repel = TRUE, cols = my_palettes[[i]]) + ggtitle(my_title))
  print(DimPlot(obj_integr, reduction = "umap", label = FALSE, group.by = interesting_idents_long[i], repel = TRUE, cols = my_palettes[[i]]) + ggtitle(my_title))
  par(mar=c(10,5,4,5))
  my_title <- paste0("Number of Cells per Cluster (", interesting_idents_short[i],")")
  print(barplot(table(obj_integr@meta.data[, interesting_idents_long[i]]), las=2, main = my_title, col = my_palettes[[i]]))
}

for(a in seq_along(levels(Idents(obj_integr)))){
  
  plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))
  mtext(levels(Idents(obj_integr))[a], side=3, line = -2, cex=3, at=-0.04, font=3,adj=0)
  seurObj <- subset(obj_integr, idents = levels(Idents(obj_integr))[a])
  for(i in seq_along(interesting_idents_long)){
    my_title <- paste0("UMAP colored by ", interesting_idents_short[i], " (",levels(Idents(obj_integr))[a],")" )
    print(DimPlot(seurObj, reduction = "umap", label = TRUE, group.by = interesting_idents_long[i],  repel = TRUE, cols = my_palettes[[i]][unique(seurObj@meta.data[, interesting_idents_long[i]])]) + ggtitle(my_title))
    print(DimPlot(seurObj, reduction = "umap", label = FALSE, group.by = interesting_idents_long[i],  repel = TRUE, cols = my_palettes[[i]][unique(seurObj@meta.data[, interesting_idents_long[i]])]) + ggtitle(my_title))
    par(mar=c(10,5,4,5))
    my_title <- paste0("Number of Cells per Cluster (", interesting_idents_short[i], ", ",levels(Idents(obj_integr))[a],")" )
    print(barplot(table(seurObj@meta.data[, interesting_idents_long[i]]), las=2, main = my_title, col = my_palettes[[i]][unique(seurObj@meta.data[, interesting_idents_long[i]])]))
  }
}
dev.off()

