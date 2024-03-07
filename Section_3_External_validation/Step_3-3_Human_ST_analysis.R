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

##########################################
### Step 3, Carcinoma region selection ###
##########################################

Add_coor <- function(data){
  data_coor <- GetTissueCoordinates(data, scale = "lowres", cols = c("imagerow", "imagecol"))
  data_coor <- P3T_coor[match(P3T@assays[["Spatial"]]@data@Dimnames[[2]], rownames(P3T_coor)),]
  data <- AddMetaData(data, data_coor$imagecol, 'col1')
  data <- AddMetaData(data, data_coor$imagerow, 'row1')
  return(data)
}

P1T <- Add_coor(P1T)
P3T <- Add_coor(P3T)
P5T <- Add_coor(P5T)
P7T <- Add_coor(P7T)
P8T <- Add_coor(P8T)
P9T <- Add_coor(P9T)
P10T <- Add_coor(P10T)
P11T <- Add_coor(P11T)

### Subset for a specific carcinoma region
P1T.sub <- subset(P1T, col1 > 130 & col1 < 320 & row1 > 290 & row1 < 470, invert = F)
P3T.sub <- subset(P3T, col1 > 240 & row1 < 400, invert = FALSE)
P5T.sub <- subset(P5T, row1 < 400 & col1 < 320)
P7T.sub <- subset(P7T, row1 > 230 & col1 < 470)
P8T.sub <- subset(P8T, col1 > 200 & col1 < 450 & row1 > 70 & row1 < 350)
P9T.sub <- subset(P9T, row1 < 300 & col1 > 290)
P10T.sub <- subset(P10T, row1 < 400 & col1 > 210)
P11T.sub <- subset(P11T, row1 > 290 & col1 > 270)

###########################################################
### Step 4, Counting numbers of RUNX2+ exhausted T cell ###
###########################################################

### Here, we counted the number of CD3D, CD3E and CD8A+ spots with 25 quantile exhausted module score
subExhT <- function(ST, idents){
  DefaultAssay(ST) <- 'Enhanced'
  a <- as.numeric(quantile(ST@assays[["Enhanced"]]@counts[1,], 0.75))
  b  <- as.numeric(quantile(ST@assays[["Enhanced"]]@counts[12,], 0.75))
  c <- as.numeric(quantile(ST@assays[["Enhanced"]]@counts[11,], 0.75))
  d <- as.numeric(quantile(ST@meta.data[["Exhaust_enhanced_module_data"]], 0.75))
  ST1 <- subset(ST, idents = idents, invert = F)
  ST1 <- subset(ST1, CD8A > a)
  ST1 <- subset(ST1, CD3D > b | CD3E > c)
  ST1 <- subset(ST1, Exhaust_enhanced_module_data > d)
  return(ST1)
}

### Number counts
P1T.subExhT <- subExhT(P1T, 'Immune/CAF') # 154
P3T.subExhT <- subExhT(P3T, 'Immune/CAF') # 64
P5T.subExhT <- subExhT(P5T, 'Immune/CAF') # 83
P7T.subExhT <- subExhT(P7T, 'Immune/CAF') # 277
P8T.subExhT <- subExhT(P8T, c('Immune','Immune/CAF')) # 446
P9T.subExhT <- subExhT(P9T, 'Immune') # 24
P10T.subExhT <- subExhT(P10T, c('Immune','Immune/CAF')) # 146
P11T.subExhT <- subExhT(P11T, 'Immune/CAF') # 428
