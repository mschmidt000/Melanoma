### reduce dimension and find clusters
### 19.05.22

# reduce dimension and plot -----------------------------------------------

map(runs_list, ~reduce_dimension_and_cluster(seurat_object = .x, figures_path = figures_path))


# save individual seurat objects ------------------------------------------

map(runs_list, ~save_seurat_objects(seurat_object = .x, output_data_path = output_data_path))
