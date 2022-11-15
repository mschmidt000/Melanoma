
import pyprojroot as here
import anndata as ad
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import scanpy as sc
import scipy
import pickle
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True
# set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

my_file = ["LE-583-KM_Analyse_Cellranger", ".loom"]
my_file = "".join(my_file)
my_path = ["data", "LE-583-KM_Analyse_Cellranger" , "velocyto", my_file]
my_path = here.here("/".join(my_path))
ldata = scv.read_loom(my_path, cache=True)

ldata_obs = pd.read_csv("cellID_Probe3.csv")
ldata_emb_cord = pd.read_csv("cell_embedding_Probe3.csv")
ldata_cell_clusters = pd.read_csv("clusters_Probe3.csv")
ldata_colors = pd.read_csv("colors_Probe3.csv")



scv.set_figure_params()


scv.pp.filter_and_normalize(adata, **params)
scv.pp.moments(adata, **params)


adata = scv.read("/home/sc.uni-leipzig.de/mf263wwia/Melanoma/data/LE-583-KM_Analyse_Cellranger/velocyto/LE-583-KM_Analyse_Cellranger.loom", cache=True)

scp -r maria@gondwanaland.izbi.uni-leipzig.de:/scratch/maria/cellranger/refdata-gex-GRCh38-2020-A ~/.
