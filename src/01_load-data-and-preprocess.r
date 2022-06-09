### load data and preprocess with seurat
### 19.05.22


# load packages and source code -------------------------------------------

pacman::p_load(Seurat, tidyverse)
source(here("src","seurat-functions.r"))


# define variables --------------------------------------------------------

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


# preprocess, filter and plot ---------------------------------------------

runs_list <- map2(runs, paths, ~create_seurat_object(
                                                data_set_name = .x,
                                                input_data_folder = .y,
                                                input_data_path = input_data_path
                                          )) %>%
                                            map(~do_qc(seurat_object = .x, figures_path = figures_path)) %>%
                                            map(~do_filtering_and_qc(seurat_object = .x, figures_path = figures_path))

# save individual seurat objects ------------------------------------------

map(runs_list, ~save_seurat_objects(seurat_object = .x, output_data_path = output_data_path))
