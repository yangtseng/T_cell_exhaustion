#################################
### Section 1, Basic analysis ###
#################################

################################################################################
### T cell extraction, cell composition and differential expression analysis ###
################################################################################

##########################################################
### T cell extraction from pre-processed seurat object ###
##########################################################

work_path = "./"
source("requirements.R")

load("murine_doubletremoval2.rds")
### It will load a pre-processed seurat object of all cells from section 1, step 1


