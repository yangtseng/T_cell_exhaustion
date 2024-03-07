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

############################
### Step 2, Module score ###
############################

### Add module score on whole ST data separately through simply average the expression level
RUNX2_module <- c('RUNX2','CTLA4','KLRK1','STAT3','IL18RAP','LGALS3','RBPJ','NRP1')
Exhaust_module <- c('PDCD1','CTLA4', 'LAG3','TIGIT','TOX', "HAVCR2")

AddAvgMeta <- function(data){
  ### Calculate the averge expression level of each module based on different slots
  Runx2_avg <- colSums(as.data.frame(data@assays$Enhanced@scale.data)[RUNX2_module, ])/8
  Exhaust_avg <- colSums(as.data.frame(data@assays$Enhanced@scale.data)[Exhaust_module, ])/6
  Runx2_avg_count <- colSums(as.data.frame(data@assays$Enhanced@counts)[RUNX2_module, ])/8
  Exhaust_avg_count <- colSums(as.data.frame(data@assays$Enhanced@counts)[Exhaust_module, ])/6
  Runx2_avg_data <- colSums(as.data.frame(data@assays$Enhanced@data)[RUNX2_module, ])/8
  Exhaust_avg_data <- colSums(as.data.frame(data@assays$Enhanced@data)[Exhaust_module, ])/6
  
  ### Setting the default assay as enhanced (result from BayesSpace)
  data <- AddMetaData(data, Runx2_avg, col.name = 'RUNX2_enhanced_module')
  data <- AddMetaData(data, Exhaust_avg, col.name = 'Exhaust_enhanced_module')
  data <- AddMetaData(data, Runx2_avg_count, col.name = 'RUNX2_enhanced_module_count')
  data <- AddMetaData(data, Exhaust_avg_count, col.name = 'Exhaust_enhanced_module_count')
  data <- AddMetaData(data, Runx2_avg_data, col.name = 'RUNX2_enhanced_module_data')
  data <- AddMetaData(data, Exhaust_avg_data, col.name = 'Exhaust_enhanced_module_data')
  
  return(data)
}

P1T <- AddAvgMeta(P1T)
P3T <- AddAvgMeta(P3T)
P5T <- AddAvgMeta(P5T)
P7T <- AddAvgMeta(P7T)
P8T <- AddAvgMeta(P8T)
P9T <- AddAvgMeta(P9T)
P10T <- AddAvgMeta(P10T)
P11T <- AddAvgMeta(P11T)

