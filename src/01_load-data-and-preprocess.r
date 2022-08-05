### load data and preprocess with seurat
### 19.05.22


# load packages and source code -------------------------------------------

pacman::p_load(Seurat, tidyverse)
source(here("src","seurat-functions.r"))
source(here("src","paths.r"))

# define variables --------------------------------------------------------

runs <- c(
  "s1_LE479_NM_T", "LE489KG_Rep", "LE493BR_Rep", "LE497NA_Rep", "LE497ST_Rep", "LE-501-DM_Rep", "LE-517-LD",
  "LE-511-MW", "LE_569_HH", "LE_577_WR", "LE_579_DE", "LE-583-KM", "LE-585-HM", "LE-589-BE_RNAseq", "LE-593-KP_RNAseq",
  "LE-595-SV_GEX", "LE-597-EG_GEX", "LE-599-SH_GEX"
  )

paths <- list.files(input_data_path, full.names = FALSE) # have to be sorted like runs!


# preprocess, filter and plot ---------------------------------------------
for(i in seq_along(runs)){
  seurat_object <- create_seurat_object(  
                      data_set_name = runs[i],
                      input_data_folder = paths[i],
                      input_data_path = input_data_path
                ) %>%
                do_qc(figures_path = figures_path)  %>%
                do_filtering_and_qc(figures_path = figures_path)
              
  save_seurat_object(seurat_object, output_data_path = output_data_path)
  rm(seurat_object)
}


