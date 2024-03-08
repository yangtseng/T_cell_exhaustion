########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 1e ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_Runx28.rds"))
