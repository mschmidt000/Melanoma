### receptor-ligand analysis for every single data set and the complete integrated data set
### 25.05.22

library(liana)
source(here("src", "liana-functions.r"))
filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)

# obj_integr@meta.data[, "predicted.id.Cell_type_nicknames_subcl_mels_nevus_marked"] <- as.character(obj_integr@meta.data[, "predicted.id.Cell_type_nicknames_subcl_mels"])
# # obj_integr@meta.data[colnames(obj_mel), interesting_idents_clust_mels] <- paste(as.character(obj_mel@meta.data[, interesting_idents_long[1]]), obj_mel@meta.data[, "integrated_snn_res.1"], sep = "_")
# nevus <- grep("LE497NA_Rep", obj_integr@meta.data[, "orig.ident"])
# obj_integr@meta.data[nevus, "predicted.id.Cell_type_nicknames_subcl_mels_nevus_marked"] <- paste(as.character(obj_integr@meta.data[nevus, "predicted.id.Cell_type_nicknames_subcl_mels"]), sep = "_", "nev")

# obj_integr@meta.data[, "predicted.id.Cell_type_nicknames_nevus_marked"] <- as.character(obj_integr@meta.data[, "predicted.id.Cell_type_nicknames"])
# nevus <- grep("LE497NA_Rep", obj_integr@meta.data[, "orig.ident"])
# obj_integr@meta.data[nevus, "predicted.id.Cell_type_nicknames_nevus_marked"] <- paste(as.character(obj_integr@meta.data[nevus, "predicted.id.Cell_type_nicknames"]), sep = "_", "nev")
# 
# obj_integr@meta.data[, "predicted.id.Cell_type_nicknames_mels_tsoi-classes"] <- as.character(obj_integr@meta.data[, "predicted.id.Cell_type_nicknames"])
# obj_integr@meta.data[, "predicted.id.Cell_type_nicknames_mels_tsoi-classes"] <- gsub("Melanocyte", "Mel", obj_integr@meta.data[, "predicted.id.Cell_type_nicknames_mels_tsoi-classes"])
# mels <- which(obj_integr@meta.data[, "predicted.id.Cell_type_nicknames"] %in% c("Melanocyte", "Mast cell"))
# # mels <- grep("Mel", obj_integr@meta.data[, "predicted.id.Cell_type_nicknames"])
# tsoi <- gsub("Tsoi_", "", obj_integr@meta.data[mels, "Tsoi_classification"])
# tsoi <- c("TrMel", "Mel", "Neu", "NeuTr", "Tr", "UnNe","Un","No" )[match(tsoi, unique(tsoi))]
# obj_integr@meta.data[mels, "predicted.id.Cell_type_nicknames_mels_tsoi-classes"] <- paste(as.character(obj_integr@meta.data[mels, "predicted.id.Cell_type_nicknames_mels_tsoi-classes"]), sep = "_", tsoi)
# 
# obj_integr@meta.data[, "predicted.id.Cell_type_nicknames_mels_tsoi-classes_nevus_marked"] <- as.character(obj_integr@meta.data[, "predicted.id.Cell_type_nicknames_mels_tsoi-classes"])
# nevus <- grep("LE497NA_Rep", obj_integr@meta.data[, "orig.ident"])
# obj_integr@meta.data[nevus, "predicted.id.Cell_type_nicknames_mels_tsoi-classes_nevus_marked"] <- paste(as.character(obj_integr@meta.data[nevus, "predicted.id.Cell_type_nicknames_mels_tsoi-classes"]), sep = "_", "nev")


Idents(obj_integr) <- "predicted.id.Cell_type_nicknames_mels_tsoi-classes_nevus_marked"

liana_obj <- liana_wrap(obj_integr, assay = "RNA", parallelize = TRUE, workers = 20)

filename <- here(output_data_path, "integrated-liana-obj-Cell_type_nicknames_nevus-marked.RData")
save(liana_obj, file = filename)
liana_obj <- liana_obj[1:4] %>%
  liana_aggregate()
save(liana_obj, file = filename)

liana_subset <- liana_obj %>%
  filter(aggregate_rank < 0.05)

filename <- here(output_data_path, "rl-results_p<0.05.csv")
write.csv2(liana_subset, file = filename,
            row.names = FALSE, quote = FALSE)


filename <- here(output_data_path, "rl-results_p<0.1.xlsx")

wb <- createWorkbook()
gc()
sheet <- createSheet(wb, sheetName = "Sheet1")
addDataFrame(liana_subset, sheet)
saveWorkbook(wb, filename)



plot_lymphs_vs_mels(liana_obj, obj_naevus_neg, "integrated")

Idents(obj_integr) <- "orig.ident"

sapply(levels(Idents(obj_integr)), function(x) {
  obj_subset <- subset(obj_integr, idents = x)
  obj_subset <- ScaleData(obj_subset, verbose = FALSE, features = rownames(obj_subset)) %>%
    RunPCA(npcs = n_dims_use, verbose = FALSE)

  liana_obj <- liana_wrap(obj_subset, assay = "RNA")
  filename <- here(output_data_path, paste0(x, "-liana-obj.RData"))
  save(liana_obj, file = filename)
  liana_obj <- liana_obj %>%
    liana_aggregate()
  save(liana_obj, file = filename)

  plot_lymphs_vs_mels(liana_obj, obj_subset, x)
})
