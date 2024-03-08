########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 1c ###
#################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

### Load data
HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_Runx28.rds"))
