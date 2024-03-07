###################################################################
### Section 3, External validation of Runx2 in exhausted T cell ###
###################################################################

##############################################################
### Validation using human spatial transcriptomics dataset ###
##############################################################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

Load(paste0(work_path, "Human_ST_bayesspace.RData"))

#############################
### Step 1, Normalization ###
#############################

### Before analyzing the enhanced ST data, we first normalized the assay
Enh_norm <- function(data){
  DefaultAssay(data) <- "Enhanced"
  data <- NormalizeData(data)
  data <- ScaleData(data)
  return(data)
}

P1T <- Enh_norm(P1T)
P3T <- Enh_norm(P3T)
P5T <- Enh_norm(P5T)
P7T <- Enh_norm(P7T)
P8T <- Enh_norm(P8T)
P9T <- Enh_norm(P9T)
P10T <- Enh_norm(P10T)
P11T <- Enh_norm(P11T)
