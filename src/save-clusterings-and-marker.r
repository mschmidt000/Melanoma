

DimPlot(obj_integr, reduction = "umap",label = TRUE,  group.by = paste0("integrated_snn_res.", c(1.2, 2,3,4,5,6)), ncol = 4) & theme(legend.position="none")

library(here)
library(Seurat)
library(ggplot2)
source(here("src","paths.r"))
source(here("src","opossom-functions.r"))

filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)

for(i in c(1.2, 2:6)){

    print(i)

    Idents(obj_integr) <- paste0("integrated_snn_res.", i)

    filename <- paste0("output/dimplot_seurat-cluster-marker_integrated-snn-res.",i,".pdf")
    pdf(filename, width = 20, height = 20)

    print(DimPlot(obj_integr, label = TRUE, repel = TRUE, raster = FALSE) +
        theme(legend.position = "none"))

    dev.off()

    # seurat_cluster_marker <- FindAllMarkers(obj_integr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    # filename <- paste0("output/seurat-cluster-marker_integrated-snn-res.",i,".csv")
    # seurat_cluster_marker <- read.csv2(filename)

    # filename <- paste0("output/seurat-cluster-marker_integrated-snn-res.",i,".csv")
    # write.csv2(seurat_cluster_marker, file = filename, row.names = FALSE, quote = FALSE)
    # filename <- paste0("output/seurat-cluster-marker_integrated-snn-res.",i,".RData")
    # save(seurat_cluster_marker, file=filename) 

    # marker_list[[paste0("integrated-snn-res.", i)]] <- seurat_cluster_marker
    # seurat_cluster_marker <- ""

}

