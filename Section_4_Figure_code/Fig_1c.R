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
Load(paste0(work_path, "Human_ST_bayesspace.RData"))
