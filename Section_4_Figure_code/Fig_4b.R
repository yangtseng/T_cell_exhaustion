########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 4b ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
PE <- read.csv(paste0(work_path, "Pathway_enrichment_Tex.csv"), header = 1)
