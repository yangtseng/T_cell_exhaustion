###########################################
### Section 1, Anlaysis of murine model ###
###########################################

#######################################
### Transcription factor prediction ###
#######################################

set.seed(1234)
work_path = "./"
source("requirements.R")

HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_TF6.rds"))
