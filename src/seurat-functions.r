#' Create seurat object, set cell name with run-prefix and calculate mitochondrial and ribosomal expression percentage
#'
#' @param data_set_name sequence of characters representing the name of the data set
#' @param input_data_folder sequence of characters representing the individual input data folder name
#' @param input_data_path character vector of general input data path
#' @return seurat object
create_seurat_object <- function(data_set_name = data_set_name, input_data_folder = input_data_folder, input_data_path = input_data_path){
  
    pacman::p_load(Seurat, tidyverse, here)
  
    seurat_object <- ReadMtx(
      mtx = here(input_data_path, input_data_folder, "filtered_feature_bc_matrix", "matrix.mtx.gz"), features = here(input_data_path, input_data_folder, "filtered_feature_bc_matrix", "features.tsv.gz"),
      cells = here(input_data_path, input_data_folder,"filtered_feature_bc_matrix", "barcodes.tsv.gz")
    )
  
    # seurat_object <- load_and_remove_ambient(seurat_object, data_set_name, input_data_folder, input_data_path)
    colnames(seurat_object) <- paste(data_set_name, colnames(seurat_object), sep = "_")
    
    seurat_object <- CreateSeuratObject(counts = as.matrix(seurat_object), min.cells = 1, min.features = 3)
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
    seurat_object[["percent.rps"]] <- PercentageFeatureSet(seurat_object, pattern = "^RPS")
    seurat_object@meta.data$orig.ident <- data_set_name
    seurat_object <- find_doublets(seurat_object)
    gc()
    seurat_object
   
}


#' Plot quality metrics as pdf
#'
#' @param seurat_object seurat object
#' @param figures_path character vector of output figures path
#' @return seurat object
do_qc <- function(seurat_object = seurat_object, figures_path = figures_path){
  
  name_run <- unique(seurat_object$orig.ident)
  pacman::p_load(Seurat, tidyverse, here)
  
  filename <- here(figures_path, paste("01", name_run, "quality-check.pdf", sep = "_"))
  pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
  
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
  number_cells <- ncol(seurat_object@assays$RNA@data)
  mtext(paste0(name_run, " (", number_cells, ")"), side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)
  print(VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rps")), main = name_run)
  p1 <- ggplot(seurat_object@meta.data) +
                  geom_point(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
                  scale_color_gradient(low = "black", high = "red")
  print(p1)
  p2 <- ggplot(seurat_object@meta.data) +
    geom_point(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.rps)) +
    scale_color_gradient(low = "black", high = "green")
  print(p2)
  print(FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" , group.by = colnames(seurat_object@meta.data)[grep("DF.class",colnames(seurat_object@meta.data))]))
  print(FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(seurat_object, feature1 = "percent.rps", feature2 = "percent.mt"))
  hist(seurat_object$nCount_RNA, col = "grey80", breaks = 100)
  hist(seurat_object$nFeature_RNA, col = "grey80", breaks = 100)
  hist(na.omit(seurat_object$percent.mt), col = "grey80", breaks = 100)
  hist(na.omit(seurat_object$percent.rps), col = "grey80", breaks = 100)
  
  plot(FeaturePlot(seurat_object,
                   reduction = "umap",
                   features = c("percent.mt", "percent.rps", "nCount_RNA", "nFeature_RNA")
  ))
  plot(DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = c("seurat_clusters")))
  plot(DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = colnames(seurat_object@meta.data)[grep("DF.class",colnames(seurat_object@meta.data))]))
  
  p3 <- VariableFeaturePlot(seurat_object)
  p4 <- LabelPoints(plot = p3, points = head(VariableFeatures(seurat_object), 10), repel = TRUE)
  
  dev.off()
  gc()
  seurat_object
  
}

#' Filter seurat object and plot filtered quality metrics as pdf
#'
#' @param seurat_object seurat object
#' @param figures_path character vector of output figures path
#' @return filtered seurat object
do_filtering_and_qc <- function(seurat_object = seurat_object, figures_path = figures_path){

  name_run <- unique(seurat_object$orig.ident)
  pacman::p_load(Seurat, tidyverse, here)
  filename <- here(figures_path, paste("01", name_run, "filtered-quality-check.pdf", sep = "_"))
  pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
  
  # seurat_object <- remove_ambient_expression(seurat_object, name_run, input_data_folder, input_data_path)
  
  filtering_cutoff_list <- list(
    "min_nCount_RNA" = quantile(seurat_object$nCount_RNA, 0.1),
    "min_nFeature_RNA" = quantile(seurat_object$nFeature_RNA, 0.1),
    "max_percent.mt" = median(seurat_object$percent.mt) + sd(seurat_object$percent.mt) * 3
  )
    
  # filtering out dead cells
  seurat_object <- subset(seurat_object,
                                 subset = nCount_RNA < filtering_cutoff_list[["min_nCount_RNA"]] & percent.mt > filtering_cutoff_list[["max_percent.mt"]] |
                                   nFeature_RNA < filtering_cutoff_list[["min_nFeature_RNA"]] & percent.mt > filtering_cutoff_list[["max_percent.mt"]] |
                                   percent.mt > filtering_cutoff_list[["max_percent.mt"]] & percent.rps < filtering_cutoff_list[["max_percent.rps"]], invert = TRUE
  )
  
  # filtering out doublets
  Idents(seurat_object) <- colnames(seurat_object@meta.data)[grep("DF.class",colnames(seurat_object@meta.data))]
  seurat_object <- subset(seurat_object, idents = "Singlet")
  Idents(seurat_object) <- "seurat_clusters"
  
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
  number_cells <- ncol(seurat_object@assays$RNA@data)
  mtext(paste0(name_run, " (", number_cells, ")"), side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)
  print(VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rps")), main = name_run)
  p1 <- ggplot(seurat_object@meta.data) +
    geom_point(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
    scale_color_gradient(low = "black", high = "red")
  print(p1)
  p2 <- ggplot(seurat_object@meta.data) +
    geom_point(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.rps)) +
    scale_color_gradient(low = "black", high = "green")
  print(p2)
  print(FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = colnames(seurat_object@meta.data)[grep("DF.class",colnames(seurat_object@meta.data))]))
  print(FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(seurat_object, feature1 = "percent.rps", feature2 = "percent.mt"))
  hist(seurat_object$nCount_RNA, col = "grey80", breaks = 100)
  hist(seurat_object$nFeature_RNA, col = "grey80", breaks = 100)
  hist(na.omit(seurat_object$percent.mt), col = "grey80", breaks = 100)
  hist(na.omit(seurat_object$percent.rps), col = "grey80", breaks = 100)
  
  plot(FeaturePlot(seurat_object,
                   reduction = "umap",
                   features = c("percent.mt", "percent.rps", "nCount_RNA", "nFeature_RNA")
  ))
  plot(DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = c("seurat_clusters")))
  plot(DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = colnames(seurat_object@meta.data)[grep("DF.class",colnames(seurat_object@meta.data))]))
  
  seurat_object <- NormalizeData(seurat_object) %>%
    FindVariableFeatures(do.plot = T, verbose = F) %>%
    ScaleData()
  p3 <- VariableFeaturePlot(seurat_object)
  p4 <- LabelPoints(plot = p3, points = head(VariableFeatures(seurat_object), 10), repel = TRUE)
  print(p4)
  
  dev.off()
  gc()
  seurat_object
    
}

#' Save seurat object
#'
#' @param seurat_object seurat object
#' @param output_data_path character vector of output data path
#' @return seurat_object
save_seurat_object <- function(seurat_object = seurat_object, output_data_path = output_data_path){
    
    name_run <- unique(seurat_object$orig.ident)
    filename <- here(output_data_path, paste0(name_run, "_seurat-obj.RData"))
    save(seurat_object, file = filename)
    gc()
    
    NULL
}

#' Load seurat object
#'
#' @param data_set_name name of to dataset to load
#' @param output_data_path character vector of output data path
#' @return seurat_object
load_seurat_object <- function(data_set_name = data_set_name, output_data_path = output_data_path){
  
  filename <- here(output_data_path, paste0(data_set_name, "_seurat-obj.RData"))
  load(filename)
  
  seurat_object
}


#' Finds doublets
#'
#' @param seurat_object seurat object
#' @return filtered seurat object
find_doublets <- function(seurat_object = seurat_object){
  
  library(DoubletFinder)
  multiplet_rate_table <- tibble(rate = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6),
                           cell_recovered_1 = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
                           cell_recovered_2 = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 1000000))
  
  multiplet_rate_table <- multiplet_rate_table %>%
                      mutate(is_between = ncol(seurat_object) >= cell_recovered_1 & ncol(seurat_object) <= cell_recovered_2)
  
  multiplet_rate <- multiplet_rate_table$rate[which(multiplet_rate_table$is_between)] * 0.01
  
  pacman::p_load(Seurat, tidyverse, here, DoubletFinder)
  
  n_dims_use <- 30
  
  perplexity <- sqrt(ncol(seurat_object@assays$RNA@counts))
  seurat_object <- seurat_object %>%
                      NormalizeData() %>%
                      FindVariableFeatures() %>%
                      ScaleData() %>%
                      RunPCA() %>%
                      RunUMAP(dims = 1:n_dims_use) %>%
                      FindNeighbors(dims = 1:n_dims_use) %>%
                      FindClusters(resolution = 1.2)
  
  sweep.res.list <- paramSweep_v3(seurat_object, PCs = 1:10, sct = FALSE)
  sweep.stat <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stat)
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  annotations <- seurat_object$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(multiplet_rate*ncol(seurat_object))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seurat_object <- doubletFinder_v3(
    seurat_object,
    PCs = 1:n_dims_use,
    pN = 0.25,
    pK = pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  ) 
  
  gc()
  
  seurat_object
  
}


#' Remove free mRNA contamination of cell 
#'
#' @param seurat_object seurat object
#' @param data_set_name sequence of characters representing the name of the data set
#' @param input_data_folder sequence of characters representing the individual input data folder name
#' @param input_data_path character vector of general input data path
#' @return filtered data matrix
load_and_remove_ambient <- function(seurat_object = seurat_object, data_set_name = data_set_name, input_data_folder = input_data_folder, input_data_path = input_data_path){
  
  pacman::p_load(Seurat, here, SoupX)
  clusters <- read_csv(here(input_data_path, input_data_folder, "graphclust", "clusters.csv"), col_types = cols(Cluster = col_character()) )
  clusters <- pull(clusters, Cluster)
  toc <- Read10X(here(input_data_path, input_data_folder, "filtered_feature_bc_matrix"))
  tod <- Read10X(here(input_data_path, input_data_folder, "raw_feature_bc_matrix"))
  sc <- SoupChannel(tod, toc)
  sc <- setClusters(sc,clusters)
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  
  gc()
  out

}

#' Dimension reduction and clustering
#'
#' @param seurat_object seurat object
#' @param figures_path character vector of output figures path
#' @return seurat object
reduce_dimension_and_cluster <- function(seurat_object = seurat_object, figures_path = figures_path){
  
  pacman::p_load(Seurat, tidyverse, here, clustree)
  
  n_dims_use <- 30
  
  perplexity <- sqrt(ncol(seurat_object))
  seurat_object <- seurat_object %>%
    RunPCA(features = VariableFeatures(seurat_object)) %>%
    FindNeighbors(dims = 1:n_dims_use) %>%
    FindClusters(resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1)) %>%
    RunTSNE(dims = 1:n_dims_use, perplexity = perplexity) %>%
    RunUMAP(dims = 1:n_dims_use) %>%
    CellCycleScoring(
      s.features = cc.genes$s.genes,
      g2m.features = cc.genes$g2m.genes,
      set.ident = FALSE
    )
  
  name_run <- unique(seurat_object$orig.ident)

  filename <- here(figures_path, paste("02", name_run, "dimension-reduction-and-clustering.pdf", sep = "_"))
  pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
  
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
  mtext(name_run, side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)
  
  plot(VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca"))
  plot(DimPlot(seurat_object, reduction = "pca"))
  plot(ElbowPlot(seurat_object))
  plot(FeaturePlot(seurat_object,
                   reduction = "umap",
                   features = c("percent.mt", "percent.rps", "nCount_RNA", "nFeature_RNA")
  ))
  plot(DimPlot(seurat_object, reduction = "tsne", label = TRUE, group.by = c("seurat_clusters", "Phase")))
  plot(DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = c("seurat_clusters", "Phase")))
  clustree(x = seurat_object@meta.data, prefix = "RNA_snn_res.") 
  plot(DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = colnames(seurat_object@meta.data)[grep("DF.class",colnames(seurat_object@meta.data))]))


  dev.off()
  
  gc()
  
  seurat_object
  
}



#' Transfer cell annotation from reference data to query data with seurat
#'
#' @param seurat_object query object
#' @param reference_object reference object
#' @param figures_path character vector of output figures path
#' @return seurat object
transfer_labels <- function(seurat_object = seurat_object, reference_object = reference_object, figures_path = figures_path){
    
    n_dims_use <- 30
    name_run <- unique(seurat_object$orig.ident)
  
    reynolds_predictions_exact <- list()
  
    my_anchors <- FindTransferAnchors(
      reference = reference_object, query = seurat_object,
      dims = 1:n_dims_use, reference.reduction = "pca"
    )
    
    predictions <- TransferData(anchorset = my_anchors, refdata = reference_object$Cell_type, dims = 1:n_dims_use)
    seurat_object$predicted.id.Cell_type <- predictions$predicted.id
    seurat_object <- AddMetaData(seurat_object, metadata = predictions)
    
    predictions <- TransferData(anchorset = my_anchors, refdata = reference_object$Cell_group, dims = 1:n_dims_use)
    seurat_object$predicted.id.Cell_group <- predictions$predicted.id
    reynolds_predictions_exact[["Cell_group"]] <- predictions
    
    predictions <- TransferData(anchorset = my_anchors, refdata = reference_object$Flow_gate, dims = 1:n_dims_use)
    seurat_object$predicted.id.Flow_gate <- predictions$predicted.id
    reynolds_predictions_exact[["Flow_gate"]] <- predictions
  
  
    interesting_idents_long <- c("predicted.id.Cell_type", "predicted.id.Cell_group", "predicted.id.Flow_gate")
    interesting_idents_short <- c("Reynolds 2021 Cell-Type", "Reynolds 2021 Cell-Group", "Reynolds 2021 Flow Gate")
  
    filename <- here(figures_path, paste("03", name_run, "cell-type-predictions.pdf", sep = "_"))
    pdf(filename, height = 7, width = 14)
    
    for (i in seq_along(interesting_idents_long)) {
      
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
      mtext(name_run, side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)
      my_title <- paste0("UMAP colored by ", interesting_idents_short[i], " (", name_run, ")")
      my_palette <- pals::polychrome(length(unique(seurat_object@meta.data[, interesting_idents_long[i]])))
      names(my_palette) <- sort(unique(seurat_object@meta.data[, interesting_idents_long[i]]))
      print(DimPlot(seurat_object, reduction = "umap", label = TRUE, cols = my_palette, group.by = interesting_idents_long[i], pt.size = 1, repel = TRUE) + ggtitle(my_title))
      par(mar = c(10, 5, 4, 5))
      my_title <- paste0("Number of Cells per Cluster (", interesting_idents_short[i], ", ", name_run, ")")
      print(barplot(table(seurat_object@meta.data[, interesting_idents_long[i]]), las = 2, main = my_title, col = my_palette[names(table(seurat_object@meta.data[, interesting_idents_long[i]]))]))
    
    }
    
    dev.off()
    
    seurat_object
  
}

#' Classify cells by Module Scores of a given Pattern
#'
#' @param seurat_object query object
#' @param name pattern which is in the name of the gene sets / module scores
#' @param set.ident setting idents by the created classification
#' @return seurat object
classify_cells <- function(
  object,
  name,
  ctrl = NULL,
  set.ident = FALSE,
  ...
) {
  if (is.null(x = ctrl)) {
    ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  }
  # object.cc <- AddModuleScore(
  #   object = object,
  #   features = features,
  #   name = name,
  #   ctrl = ctrl,
  #   search = TRUE
  # )
  cc.columns <- grep(pattern = name, x = colnames(x = object[[]]), value = TRUE)
  cc.scores <- object[[cc.columns]]
  CheckGC()
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores, names_scores = names(scores)) {
      if (all(scores < 0)) {
        return("No")
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undecided')
        } else {
          return(gsub('[[:digit:]]+', '',names_scores[which(x = scores == max(scores))]))
        }
      }
    }
  )
  cc.scores_new <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores_new) <- c('rownames', colnames(cc.scores), paste(name, "classification", sep = "_"))
  rownames(x = cc.scores_new) <- cc.scores_new$rownames
  cc.scores_new <- cc.scores_new[, -1]
  object[[colnames(cc.scores_new)]] <- cc.scores_new
  if (set.ident) {
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- paste(name, "classification", sep = "_")
  }
  return(object)
}

#' Annotate seurat clusters by reference clusters
#'
#' @param seurat_object query object
#' @param reference_clusters reference clusters by which the seurat clusters will be annotated
#' @return seurat object
annotate_clusters <- function(
  object,
  reference_clusters
) {
  tab <- table(object$seurat_clusters, object@meta.data[,reference_clusters])
  CheckGC()
  annotated_cluster <- apply(
    X = tab,
    MARGIN = 1,
    FUN = function(cluster) {
      if (all(cluster == 0)) {
        return("No")
      } else {
          return(names(cluster)[which(x = cluster == max(cluster))])
      }
    }
  )
  return(annotated_cluster)
}
