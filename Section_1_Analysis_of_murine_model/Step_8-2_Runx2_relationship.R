###########################################
### Section 1, Anlaysis of murine model ###
###########################################

############################################################
### Runx2 relationship with markers of T cell exhaustion ###
############################################################

set.seed(1234)
work_path = "./"
source("requirements.R")

HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_Pathway7.rds"))
