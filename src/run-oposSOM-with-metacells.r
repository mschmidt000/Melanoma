sters-and-tso-classes.r### run opossom with metacells
### 12.10.22
library(Seurat)
library(oposSOM)
library(here)
library(reticulate)
source(here("src","paths.r"))
source(here("src","seurat-functions.r"))
source(here("src","opossom-functions.r"))
source(here("src","convert-huge-sparse-matrix-to-matrix.r"))
# py_install("python-igraph")
# py_install("leidenalg")
# py_install("numpy")

filename <- here(output_data_path, "integrated-seurat-obj.RData")
load(filename)

annotated_clusters <- annotate_clusters(obj_integr, "predicted.id.Cell_type_nicknames")
obj_integr$seurat_clusters_annotated <- annotated_clusters[obj_integr$seurat_clusters]
obj_integr$seurat_clusters_annotated <- paste0(obj_integr$seurat_clusters_annotated,"_c",obj_integr$seurat_clusters )
obj_integr$seurat_clusters_annotated <- factor(obj_integr$seurat_clusters_annotated, levels = unique(obj_integr$seurat_clusters_annotated)[order(unique(obj_integr$seurat_clusters))])

labels <- paste0(obj_integr$seurat_clusters_annotated, "_", obj_integr@meta.data$orig.ident)
names(labels) <- colnames(obj_integr)
labels <- gsub(" ", "",labels)

# ## remove meta-cells with less than 5 cells ###

# tab <- table(labels)
# labels <- labels[which(labels %in% names(table(labels)[table(labels)>=5]))]
# barplot(tab, las=2, log="y")

# ## subclustering-determine number of subclusters ###

# labels.clusterNo <- ceiling(sort(table(labels)) / 20)
# barplot(labels.clusterNo, las = 2)

# labels.clusterNo.interest <- names(labels.clusterNo)[grep("Melanocyte",names(labels.clusterNo))]
# labels.clusterNo[which(labels.clusterNo>2)] <- 2
# labels.clusterNo[grep("Melanocyte", names(labels.clusterNo))] <- labels.clusterNo.interest

# barplot(labels.clusterNo, las=2)

# o <- order(sapply(strsplit(names(labels.clusterNo), "_"), "[", 1))
# labels.clusterNo <- labels.clusterNo[o]
# o <- order(sapply(strsplit(names(labels.clusterNo), "_"), function(x){ as.numeric(substr(x[2], 2, nchar(x[2])))}))
# labels.clusterNo <- labels.clusterNo[o]

# barplot(labels.clusterNo, las=2)

# ## subclustering part 1 ###

# metacell.labels <- rep(NA, ncol(obj_integr))
# names(metacell.labels) <- colnames(obj_integr)
# metacell.labels[names(labels)] <- labels

# obj_integr[["metacellLabels"]] <- metacell.labels
# obj_integr[["cellInMetacell"]] <- !is.na(metacell.labels)

# metacell.data <- matrix(NA, nrow(obj_integr@assays$RNA@data), 0, dimnames=list(rownames(obj_integr@assays$RNA@data),c()))

# ## subclustering part 2 ###

# pb <- txtProgressBar(max = length(labels.clusterNo), style = 3)
# for( x in names(labels.clusterNo)){

#   mc.cells <- names(which(metacell.labels==x))
#   expr <- data.matrix(obj_integr@assays$RNA@data[,mc.cells])

#   km <- kmeans(t(expr), centers = labels.clusterNo[x])
#   lab <- paste0(x, "_k", seq(max(km$cluster)))

#   expr <- t(km$centers)
#   colnames(expr) <- lab

#   metacell.data <- cbind(metacell.data, expr)

#   obj_integr$metacellLabels[mc.cells] <- paste0(obj_integr$metacellLabels[mc.cells], " x", km$cluster)

#   setTxtProgressBar(pb, pb$getVal()+1)

# }
# pf$kill()

load("output/metacell-indata-opossom.RData")

env <- opossom.new(list(
  dataset.name = paste(Sys.Date(), "melanoma", "metacells", sep = "_"),
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
env$indata <- metacell.data
env$group.labels <- sapply(strsplit(colnames(metacell.data), "_"), function(x){ paste(x[1], x[2], sep = "_") })
names(env$group.labels) <- colnames(metacell.data)

env$group.colors <- env$group.labels %>%
  dplyr::n_distinct() %>%
  randomcoloR::distinctColorPalette()
names(env$group.colors) <- unique(env$group.labels) 
env$group.colors <- env$group.colors[match(env$group.labels, names(env$group.colors))]
names(env$group.colors) <- colnames(metacell.data)
barplot(1:length(env$group.colors)*2, col = env$group.colors, border = NA)

env$group.labels <- env$group.labels[sort(names(env$group.labels))]
str(env$group.labels)
env$group.colors <- env$group.colors[names(env$group.labels)]
str(env$group.colors)
barplot(1:length(env$group.colors)*2, col = env$group.colors, border = NA)

env$indata <- env$indata[, names(env$group.labels)]
str(env$indata)

opossom.run(env)