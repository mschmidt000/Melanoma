#' Create seurat object, set cell names with run-prefix and calculate mitochondrial and ribosomal expression percentage
#'
#' @param data_set_name sequence of characters representing the name of the data set
#' @param input_data_folder sequence of characters representing the individual input data folder name
#' @param input_data_path character vector of general input data path
#' @return seurat object
create_seurat_object <- function(data_set_name = data_set_name, input_data_folder = input_data_folder, input_data_path = input_data_path){
    pacman::p_load(Seurat, tidyverse, here)
    seurat_object <- ReadMtx(
      mtx = here(input_data_path, input_data_folder, "matrix.mtx.gz"), features = here(input_data_path, input_data_folder, "features.tsv.gz"),
      cells = here(input_data_path, input_data_folder, "barcodes.tsv.gz")
    )
    colnames(seurat_object) <- paste(data_set_name, colnames(seurat_object), sep = "_")
    seurat_object <- CreateSeuratObject(counts = as.matrix(seurat_object), min.cells = 1, min.features = 1)
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
    seurat_object[["percent.rps"]] <- PercentageFeatureSet(seurat_object, pattern = "^RPS")
    seurat_object@meta.data$orig.ident <- data_set_name

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
    print(FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt"))
    print(FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "percent.mt"))
    print(FeatureScatter(seurat_object, feature1 = "percent.rps", feature2 = "percent.mt"))
    hist(seurat_object$nCount_RNA, col = "grey80", breaks = 100)
    hist(seurat_object$nFeature_RNA, col = "grey80", breaks = 100)
    hist(na.omit(seurat_object$percent.mt), col = "grey80", breaks = 100)
    hist(na.omit(seurat_object$percent.rps), col = "grey80", breaks = 100)
 dev.off()
 
 seurat_object
  
}

#' Filter seurat objects and plot filtered quality metrics as pdf
#'
#' @param seurat_object seurat object
#' @param figures_path character vector of output figures path
#' @return filtered seurat object
do_filtering_and_qc <- function(seurat_object = seurat_object, figures_path = figures_path){

  name_run <- unique(seurat_object$orig.ident)
  pacman::p_load(Seurat, tidyverse, here)
  filename <- here(figures_path, paste("01", name_run, "filtered-quality-check.pdf", sep = "_"))
  pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
  
  filtering_cutoff_list <- list(
    "min_nCount_RNA" = quantile(seurat_object$nCount_RNA, 0.25),
    "max_nCount_RNA" = median(seurat_object$nCount_RNA) + sd(seurat_object$nCount_RNA) * 3,
    "max_percent.mt" = median(seurat_object$percent.mt) + sd(seurat_object$percent.mt) * 3,
    "max_percent.rps" = quantile(seurat_object$percent.rps, 0.30),
    "min_nFeature_RNA" = quantile(seurat_object$nFeature_RNA, 0.25)
  )
    
  # filtering out dead cells
  seurat_object <- subset(seurat_object,
                                 subset = nCount_RNA < filtering_cutoff_list[["min_nCount_RNA"]] & percent.mt > filtering_cutoff_list[["max_percent.mt"]] |
                                   nFeature_RNA < filtering_cutoff_list[["min_nFeature_RNA"]] & percent.mt > filtering_cutoff_list[["max_percent.mt"]] |
                                   percent.mt > filtering_cutoff_list[["max_percent.mt"]] & percent.rps < filtering_cutoff_list[["max_percent.rps"]], invert = TRUE
  )
  
  # filtering out doublets
  seurat_object <- subset(seurat_object,
                                 subset = nCount_RNA < filtering_cutoff_list["max_nCount_RNA"]
  )
  
  
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
  number_cells <- ncol(seurat_object@assays$RNA@data)
  mtext(paste0(name_run, " (", number_cells, ")"), side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)
  print(VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")), main = name_run)
  print(FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(seurat_object, feature1 = "percent.rps", feature2 = "percent.mt"))
  hist(seurat_object$nCount_RNA, col = "grey80", breaks = 100)
  hist(seurat_object$nFeature_RNA, col = "grey80", breaks = 100)
  hist(na.omit(seurat_object$percent.mt), col = "grey80", breaks = 100)
  hist(na.omit(seurat_object$percent.rps), col = "grey80", breaks = 100)

  seurat_object <- NormalizeData(seurat_object) %>%
    FindVariableFeatures(do.plot = T, verbose = F) %>%
    ScaleData()
  plot1 <- VariableFeaturePlot(seurat_object)
  plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(seurat_object), 10), repel = TRUE)
  pplot2 + ggtitle(name_run)
  
  dev.off()
  
  seurat_object
    
}

#' Save seurat objects in list
#'
#' @param seurat_object seurat object
#' @param output_data_path character vector of output data path
#' @return seurat_object
save_seurat_objects <- function(seurat_object = seurat_object, output_data_path = output_data_path){
    
    name_run <- unique(seurat_object$orig.ident)
    filename <- here(output_data_path, paste0(name_run, "_seurat-obj.RData"))
    save(seurat_object, file = filename)

}
