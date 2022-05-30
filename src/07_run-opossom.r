### run opossom on melanocyte cells
### 26.05.22
library(oposSOM)

filename <- here(output_data_path, "integrated-seurat-obj-melanocytes.RData")
load(filename)
Idents(obj_mel) <- "integrated_snn_res.1"
obj_mel_ds <- subset(obj_mel, downsample = 300) %>%
  ScaleData(verbose = FALSE, features = rownames(obj_mel))
rm(obj_mel)
env <- opossom.new(list(
  dataset.name = "2022-05-26-melanocytes",
  dim.1stLvlSom = "auto",
  dim.2ndLvlSom = "auto",
  activated.modules = list(
    "reporting" = TRUE,
    "primary.analysis" = TRUE,
    "sample.similarity.analysis" = FALSE,
    "geneset.analysis" = TRUE,
    "psf.analysis" = TRUE,
    "group.analysis" = TRUE,
    "difference.analysis" = TRUE
  ),
  standard.spot.modules = "group.overexpression"
))


# definition of indata, group.labels and group.colors
env$indata <- GetAssayData(obj_mel_ds, slot = "scale.data", assay = "integrated") %>%
  as.matrix()
env$group.labels <- paste("Melanocyte", obj_mel_ds$integrated_snn_res.1, sep = "_")
env$group.colors <- obj_mel_ds$integrated_snn_res.1 %>%
  n_distinct() %>%
  pals::glasbey()
ord <- order(env$group.labels)
env$indata <- env$indata[, ord]
env$group.labels <- env$group.labels[ord]
env$group.colors <- env$group.colors[ord]
rm(obj_mel_ds)
opossom.run(env)
