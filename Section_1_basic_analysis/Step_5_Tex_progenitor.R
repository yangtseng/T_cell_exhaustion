#################################
### Section 1, Basic analysis ###
#################################

#######################################################
### Further analysis of Tex progenitor in cluster 4 ###
#######################################################

set.seed(1234)
work_path = "./"
source("requirements.R")

load("murine_tcell3.rds")
### It will load a pre-processed seurat object of T cells from section 1, step 3

################################
### Subset of Tex progenitor ###
################################
