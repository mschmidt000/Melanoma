
library(dyplr)
library(ggplot2)
DefaultAssay(obj_integr) <- "integrated"
obj_integr <- FindClusters(obj_integr, resolution = c(2,3,4,5,6))

DimPlot(obj_integr, reduction = "umap", group.by = paste0("integrated_snn_res.", c(1.2, 2,3,4,5,6)), ncol = 4) & theme(legend.position="none")

clustree(x = obj_integr@meta.data, prefix = "integrated_snn_res.")

DefaultAssay(obj_integr) <- "RNA"

marker_list <- list()
pdf("output/UMAP-of-all-clusterings.pdf")
for(i in c(1.2, 2:6)){

    print(i)

    Idents(obj_integr) <- paste0("integrated_snn_res.", i)

    DimPlot(obj_integr, label = TRUE, repel = TRUE, raster = FALSE) +
        theme(legend.position = "none") 

    seurat_cluster_marker <- FindAllMarkers(obj_integr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    filename <- paste0("seurat-cluster-marker_integrated-snn-res.",i,".csv")
    write.csv2(seurat_cluster_marker, file = filename, row.names = FALSE, quote = FALSE)
    filename <- paste0("seurat-cluster-marker_integrated-snn-res.",i,".RData")
    save(seurat_cluster_marker, file="BF-01-02-03-all-Analysis_seuratObj_integrated_seurat.cluster.markers_percent.mt<15.RData") 

    marker_list[[paste0("integrated-snn-res.", i)]] <- seurat_cluster_marker
    seurat_cluster_marker <- ""

}
dev.off()