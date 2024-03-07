###################################################################
### Section 3, External validation of Runx2 in exhausted T cell ###
###################################################################

##############################################################
### Validation using human spatial transcriptomics dataset ###
##############################################################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

##############################
### Step 1, Basic analysis ###
##############################

### Load ST data
P1T <- readRDS(paste0(work_path, "Human_ST/P1T_Spatial.rds.gz"))
P3T <- readRDS(paste0(work_path, "Human_ST/P3T_Spatial.rds.gz"))
P5T <- readRDS(paste0(work_path, "Human_ST/P5T_Spatial.rds.gz"))
P7T <- readRDS(paste0(work_path, "Human_ST/P7T_Spatial.rds.gz"))
P8T <- readRDS(paste0(work_path, "Human_ST/P8T_Spatial.rds.gz"))
P9T <- readRDS(paste0(work_path, "Human_ST/P9T_Spatial.rds.gz"))
P10T <- readRDS(paste0(work_path, "Human_ST/P10T_Spatial.rds.gz"))
P11T <- readRDS(paste0(work_path, "Human_ST/P11T_Spatial.rds.gz"))
