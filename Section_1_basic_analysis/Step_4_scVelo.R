###########################################
### Section 1, Anlaysis of murine model ###
###########################################

###############################
### RNA velocity via scVelo ###
###############################

set.seed(1234)
work_path = "./"
source("requirements.R")

load("murine_tcell_modulescore4.rds")

### The scVelo is running under python 
### Here, we converted original Seurat object into h5ad for scanpy and scVelo in python

###################################
### Step 1, write H5Seurat file ###
###################################

SaveH5Seurat(HCC.tcell, paste0(workt_path, "HCC_tcell.H5Seurat"), overwrite = TRUE)

###############################################
### Step 2, convert H5Seurat file into h5ad ###
###############################################

Convert(paste0(work_path, "HCC_tcell.h5Seurat"), dest = "h5ad")

### Please refer to Step_4_scVelo.py for later analysis
