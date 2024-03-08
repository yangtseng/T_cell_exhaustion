########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 2b ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds"))
