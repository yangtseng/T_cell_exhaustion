########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 5d ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
load("./external_validation/Human_ST_subT.RData")
### Including P1T.subT, P3T.subT, P5T.subT, P7T.subT, P8T.subT, P9T.subT, P10T.subT and P11T.subT
