###################################################
### Section 2, Anlaysis of cell line like model ###
###################################################

#######################################
### Transcription factor prediction ###
#######################################

set.seed(1234)
work_path = "./"
source("requirements.R")

tcell <- readRDS(paste0(work_path, "tcell_multiome2.rds"))

### We mainly followed the tutorial of SCENIC 
### http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

#############################################
### Step 1, Initialize settings of SCENIC ###
#############################################

org <- "mgi" # or hgnc, or dmel
dbDir <- paste0(work_path, "cisTarget-databases") # RcisTarget databases location
myDatasetTitle <- "TF_in_cellline_like_model" # choose a name for your analysis
