###################################################################
### Section 3, External validation of Runx2 in exhausted T cell ###
###################################################################

##############################################################
### Validation using human spatial transcriptomics dataset ###
##############################################################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

###################################
### Step 1, Spatial enhancement ###
###################################

### Load ST data
load(paste(work_path, "Human_ST_annotated.RData"))
