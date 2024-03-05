###########################################
### Section 1, Anlaysis of murine model ###
###########################################

#######################################
### Transcription factor prediction ###
#######################################

set.seed(1234)
work_path = "./"
source("requirements.R")

### We mainly followed the tutorial of SCENIC 
### http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

######################################
### Step 1, Import GRNboost output ###
######################################

GRNBoost_output <- read.delim(paste0(work_path, "int/1.1_grn_output.tsv"), header=FALSE)
colnames(GRNBoost_output) <- c("TF","Target","weight")
saveRDS(GRNBoost_output, file="int/1.4_GENIE3_linkList.Rds")
