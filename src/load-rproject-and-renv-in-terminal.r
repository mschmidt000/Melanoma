install.packages("renv")

# load the renv library + project
library(renv)
setwd("Melanoma")
renv::load()
renv::restore()
