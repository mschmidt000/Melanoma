### predict cell types with data set from reynolds 2021
### 20.05.22

cell_type_to_exclude <- "Schwann"

filename <- here(literature_path, "Reynolds2021_Data-all_seuratObj.RData")
load(filename)
Idents(seurObj) <- "Cell_type"
seurObj <- subset(seurObj, idents = cell_type_to_exclude, invert = TRUE)

reynolds_predictions_exact <- list()

for(i in seq_along(runs_list)){
  
  reynolds_predictions_exact[[runs[i]]] <- list()
  
  my_anchors <- FindTransferAnchors(reference = seurObj, query = runs_list[[runs[i]]],
                                    dims = 1:n_dims_use, reference.reduction = "pca")
  
  predictions <- TransferData(anchorset = my_anchors, refdata = seurObj$Cell_type, dims = 1:n_dims_use)
  runs_list[[runs[i]]]$predicted.id.Cell_type <- predictions$predicted.id
  reynolds_predictions_exact[[runs[i]]][["Cell_type"]] <- predictions
  
  predictions <- TransferData(anchorset = my_anchors, refdata = seurObj$Cell_group, dims = 1:n_dims_use)
  runs_list[[runs[i]]]$predicted.id.Cell_group <- predictions$predicted.id
  reynolds_predictions_exact[[runs[i]]][["Cell_group"]] <- predictions
  
  predictions <- TransferData(anchorset = my_anchors, refdata = seurObj$Flow_gate, dims = 1:n_dims_use)
  runs_list[[runs[i]]]$predicted.id.Flow_gate <- predictions$predicted.id
  reynolds_predictions_exact[[runs[i]]][["Flow_gate"]] <- predictions

}


save(reynolds_predictions_exact, file = here(output_data_path,"reynolds-predictions-exact.RData"))
rm(reynolds_predictions_exact)
rm(seurObj)

for(i in seq_along(runs_list)){
  seurat_obj <- runs_list[[runs[i]]]
  filename <- here(output_data_path,paste0(runs[i], "_seurat-obj.RData"))
  save(seurat_obj, file = filename)
}
rm(seurat_obj)


interesting_idents_long <- c("predicted.id.Cell_type", "predicted.id.Cell_group", "predicted.id.Flow_gate" )
interesting_idents_short <- c("Reynolds 2021 Cell-Type", "Reynolds 2021 Cell-Group", "Reynolds 2021 Flow Gate")


filename <- here(figures_path, "03-cell-type-predictions.pdf")
pdf(filename, height = 7, width = 14)

for(a in seq_along(runs_list)){
  
  seurat_obj <- runs_list[[runs[a]]]
  for(i in seq_along(interesting_idents_long)){
    
      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))
      mtext(names(runs_list)[a], side=3, line = -2, cex=3, at=-0.04, font=3,adj=0)
      my_title <- paste0("UMAP colored by ", interesting_idents_short[i], " (",names(runs_list)[a],")" )
      my_palette <- pals::polychrome(length(unique(seurat_obj@meta.data[, interesting_idents_long[i]])))
      names(my_palette) <- sort(unique(seurat_obj@meta.data[, interesting_idents_long[i]]))
      print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, cols = my_palette, group.by = interesting_idents_long[i], pt.size = 1, repel = TRUE) + ggtitle(my_title))
      par(mar=c(10,5,4,5))
      my_title <- paste0("Number of Cells per Cluster (", interesting_idents_short[i], ", ",names(runs_list)[a],")" )
      print(barplot(table(seurat_obj@meta.data[, interesting_idents_long[i]]), las=2, main = my_title, col = my_palette[names(table(seurat_obj@meta.data[, interesting_idents_long[i]]))]))
  
    }
}

dev.off()
