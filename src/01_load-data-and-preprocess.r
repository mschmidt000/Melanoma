### load data and preprocess with seurat
### 19.05.22

pacman::p_load(Seurat, tidyverse)

runs <- c(
  "s1_LE479_NM_T", "LE489KG_Rep", "LE493BR_Rep", "LE497NA_Rep", "LE497ST_Rep", "LE-501-DM_Rep", "LE-517-LD",
  "LE-511-MW", "LE_569_HH", "LE_577_WR", "LE_579_DE", "LE-583-KM", "LE-585-HM", "LE-597-EG_GEX", "LE-595-SV_GEX",
  "LE-593-KP_RNAseq", "LE-589-BE_RNAseq"
)
input_data_path <- file.path("data")
output_data_path <- file.path("output")
figures_path <- file.path("figs")
literature_path <- "literature"

paths <- list.dirs(input_data_path, full.names = FALSE)[-which(list.dirs(input_data_path, full.names = FALSE) %in% "")]


runs_list <- list()
filtering_cutoff_list <- list()

for (i in seq_along(runs)) {
  runs_list[[runs[i]]] <- list()
  runs_list[[runs[i]]] <- ReadMtx(
    mtx = here(input_data_path, paths[i], "matrix.mtx.gz"), features = here(input_data_path, paths[i], "features.tsv.gz"),
    cells = here(input_data_path, paths[i], "barcodes.tsv.gz")
  )
  colnames(runs_list[[runs[i]]]) <- paste(runs[i], colnames(runs_list[[runs[i]]]), sep = "_")
  runs_list[[runs[i]]] <- CreateSeuratObject(counts = as.matrix(runs_list[[runs[i]]]), min.cells = 1, min.features = 1)
  runs_list[[runs[i]]][["percent.mt"]] <- PercentageFeatureSet(runs_list[[runs[i]]], pattern = "^MT-")
  runs_list[[runs[i]]][["percent.rps"]] <- PercentageFeatureSet(runs_list[[runs[i]]], pattern = "^RPS")
  runs_list[[runs[i]]]@meta.data$orig.ident <- runs[i]
}

filename <- here(figures_path, "01_quality-check.pdf")
pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
for (i in seq_along(runs_list)) {
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
  number_cells <- ncol(runs_list[[runs[i]]]@assays$RNA@data)
  mtext(paste0(runs[i], " (", number_cells, ")"), side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)
  print(VlnPlot(object = runs_list[[runs[i]]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rps")), main = runs[i])
  print(FeatureScatter(runs_list[[runs[i]]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(runs_list[[runs[i]]], feature1 = "nFeature_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(runs_list[[runs[i]]], feature1 = "percent.rps", feature2 = "percent.mt"))
  hist(runs_list[[runs[i]]]$nCount_RNA, col = "grey80", breaks = 100)
  hist(runs_list[[runs[i]]]$nFeature_RNA, col = "grey80", breaks = 100)
  hist(na.omit(runs_list[[runs[i]]]$percent.mt), col = "grey80", breaks = 100)
  hist(na.omit(runs_list[[runs[i]]]$percent.rps), col = "grey80", breaks = 100)
}
dev.off()

filename <- here(figures_path, "01_filtered-quality-check.pdf")
pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
for (i in seq_along(runs_list)) {
  filtering_cutoff_list[[runs[i]]] <- c(
    "min_nCount_RNA" = quantile(runs_list[[runs[i]]]$nCount_RNA, 0.25),
    "max_nCount_RNA" = median(runs_list[[runs[i]]]$nCount_RNA) + sd(runs_list[[runs[i]]]$nCount_RNA) * 3,
    "max_percent.mt" = median(runs_list[[runs[i]]]$percent.mt) + sd(runs_list[[runs[i]]]$percent.mt) * 3,
    "max_percent.rps" = quantile(runs_list[[runs[i]]]$percent.rps, 0.30),
    "min_nFeature_RNA" = quantile(runs_list[[runs[i]]]$nFeature_RNA, 0.25)
  )

  # filtering out dead cells
  runs_list[[runs[i]]] <- subset(runs_list[[runs[i]]],
    subset = nCount_RNA < filtering_cutoff_list[[i]]["min_nCount_RNA"] & percent.mt > filtering_cutoff_list[[i]]["max_percent.mt"] |
      nFeature_RNA < filtering_cutoff_list[[i]]["min_nFeature_RNA"] & percent.mt > filtering_cutoff_list[[i]]["max_percent.mt"] |
      percent.mt > filtering_cutoff_list[[i]]["max_percent.mt"] & percent.rps < filtering_cutoff_list[[i]]["max_percent.rps"], invert = TRUE
  )

  # filtering out doublets
  runs_list[[runs[i]]] <- subset(runs_list[[runs[i]]],
    subset = nCount_RNA < filtering_cutoff_list[[i]]["max_nCount_RNA"]
  )


  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1))
  number_cells <- ncol(runs_list[[runs[i]]]@assays$RNA@data)
  mtext(paste0(runs[i], " (", number_cells, ")"), side = 3, line = -2, cex = 3, at = -0.04, font = 3, adj = 0)
  print(VlnPlot(object = runs_list[[runs[i]]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt")), main = runs[i])
  print(FeatureScatter(runs_list[[runs[i]]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(runs_list[[runs[i]]], feature1 = "nFeature_RNA", feature2 = "percent.mt"))
  print(FeatureScatter(runs_list[[runs[i]]], feature1 = "percent.rps", feature2 = "percent.mt"))
  hist(runs_list[[runs[i]]]$nCount_RNA, col = "grey80", breaks = 100)
  hist(runs_list[[runs[i]]]$nFeature_RNA, col = "grey80", breaks = 100)
  hist(na.omit(runs_list[[runs[i]]]$percent.mt), col = "grey80", breaks = 100)
  hist(na.omit(runs_list[[runs[i]]]$percent.rps), col = "grey80", breaks = 100)
}
dev.off()

filename <- here(figures_path, "01_top-10-variable-features.pdf")
pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
for (i in seq_along(runs_list)) {
  runs_list[[runs[i]]] <- NormalizeData(runs_list[[runs[i]]]) %>%
    FindVariableFeatures(do.plot = T, verbose = F) %>%
    ScaleData()
  plot1 <- VariableFeaturePlot(runs_list[[runs[i]]])
  plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(runs_list[[runs[i]]]), 10), repel = TRUE)
  print(plot1 + plot2 + ggtitle(runs[i]))
}
dev.off()

for (i in seq_along(runs_list)) {
  seurat_obj <- runs_list[[runs[i]]]
  filename <- here(output_data_path, paste0(runs[i], "_seurat-obj.RData"))
  save(seurat_obj, file = filename)
}
rm(seurat_obj)
