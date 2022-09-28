conda create --name r4.2.1-base
conda activate r4.2.1-base
conda install -c conda-forge r-base=4.2.1

# save package file from .cache of the old repo to the .cache of the new repo
cp -r /home/sc.uni-leipzig.de/mf263wwia/.cache/R/renv/cache/v5/R-4.2/x86_64-pc-linux-gnu/* /home/sc.uni-leipzig.de/mf263wwia/.cache/R/renv/cache/v5/R-4.2/x86_64-conda-linux-gnu/

cp -r ~/Melanoma/renv/library/R-4.2/x86_64-pc-linux-gnu/* ~/Melanoma/renv/library/R-4.2/x86_64-conda-linux-gnu/

