###########################################
### Section 1, Anlaysis of murine model ###
###########################################

#######################################
### Transcription factor prediction ###
#######################################

set.seed(1234)
work_path = "./Cellline/"
source("requirements.R")

### We mainly followed the tutorial of SCENIC 
### http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

### load SCENIC option
scenicOptions <- readRDS(paste0(work_path, "int/scenicOptions.rds"))
exprMat_filtered_log <- readRDS(paste0(work_path, "int/exprMat_filtered_log.rds"))

