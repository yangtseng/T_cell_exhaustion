###################################################################
### Section 3, External validation of Runx2 in exhausted T cell ###
###################################################################

##############################################################
### Validation using human spatial transcriptomics dataset ###
##############################################################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

### Load ST data
P1T <- readRDS(paste0(work_path, "Human_ST/P1T_Spatial.rds.gz"))
P3T <- readRDS(paste0(work_path, "Human_ST/P3T_Spatial.rds.gz"))
P5T <- readRDS(paste0(work_path, "Human_ST/P5T_Spatial.rds.gz"))
P7T <- readRDS(paste0(work_path, "Human_ST/P7T_Spatial.rds.gz"))
P8T <- readRDS(paste0(work_path, "Human_ST/P8T_Spatial.rds.gz"))
P9T <- readRDS(paste0(work_path, "Human_ST/P9T_Spatial.rds.gz"))
P10T <- readRDS(paste0(work_path, "Human_ST/P10T_Spatial.rds.gz"))
P11T <- readRDS(paste0(work_path, "Human_ST/P11T_Spatial.rds.gz"))

##################################
### Step 1, Data preprocessing ###
##################################

SCT_preprocessing <- function(data, dim, verbose){
  data <- SCTransform(data, assay = 'Spatial', verbose = verbose, min_cells = 3)
  data <- RunPCA(data, assay = 'SCT', verbose = verbose)
  data <- FindNeighbors(data, reduction = 'pca', dims = 1:dim)
  data <- FindClusters(data, verbose = verbose, resolution = 1.5)
  data <- RunUMAP(data, reduction = 'pca', dims = 1:dim)
  
  return(data)  
}

### Here we preprocessed the data with dims = 20 and resolution = 1.5 separately
P1T <- SCT_preprocessing(P1T, 20, F)
P3T <- SCT_preprocessing(P3T, 20, F)
P5T <- SCT_preprocessing(P5T, 20, F)
P7T <- SCT_preprocessing(P7T, 20, F)
P8T <- SCT_preprocessing(P8T, 20, F)
P9T <- SCT_preprocessing(P9T, 20, F)
P10T <- SCT_preprocessing(P10T, 20, F)
P11T <- SCT_preprocessing(P11T, 20, F)

####################################
### Step 2, Cell type annotation ###
####################################

### We manually annotated the cell type of each cluster based on the specific cell-type markers
### cell type markers: 'HAMP','SPINK1','CD3D','BANK1','ITGAX','COL1A2','CLDN5','CD80','CHIT1','CYP1A2','TOP2A','CD3E','CD19','CD33','FAP','PECAM1','CD86','VSIG4','APOA5','APOE','CD3G','CD22','CD68','DCN','CD34','CD1E','CLEC5A'
### The cell type annotation result of each sample
P1T.celltype <- c('Malignant hepatocyte','Malignant hepatocyte','Immune/CAF','Hepatocyte','Endothelial/CAF','Malignant hepatocyte','Malignant hepatocyte',
                  'Malignant hepatocyte','Malignant hepatocyte','Hepatocyte','Hepatocyte','Malignant hepatocyte','Endothelial','Malignant hepatocyte')
P3T.celltype <- c('Malignant hepatocyte','CAF','Malignant hepatocyte','Immune/CAF','Malignant hepatocyte','Endothelial/CAF','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte',
                  'Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Immune/CAF','Hepatocyte')
P5T.celltype <- c('Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Immune/CAF','Malignant hepatocyte','Malignant hepatocyte','Hepatocyte',
                  'Malignant hepatocyte','Hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Hepatocyte','Endothelial/CAF','Malignant hepatocyte','Immune/CAF')
P7T.celltype <- c('Malignant hepatocyte','Malignant hepatocyte','Immune/CAF','Immune/CAF','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Immune/CAF',
                  'Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Endothelial/CAF','Malignant hepatocyte','Endothelial/CAF')
P8T.celltype <- c('Immune/CAF','Malignant hepatocyte','Malignant hepatocyte','Immune/CAF','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','CAF','Malignant hepatocyte',
                  'Malignant hepatocyte','Endothelial/CAF','Immune','Endothelial/CAF','Immune/CAF','Malignant hepatocyte','Malignant hepatocyte','Immune','CAF','Malignant hepatocyte','CAF')
P9T.celltype <- c('Hepatocyte','Hepatocyte','Hepatocyte','Hepatocyte','Immune','Endothelial/CAF','Malignant hepatocyte','Endothelial/CAF','Malignant hepatocyte','Hepatocyte','Malignant hepatocyte','Malignant hepatocyte',
                  'Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Hepatocyte','Immune','Immune')
P10T.celltype <-  c('Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Endothelial/CAF','Malignant hepatocyte','Malignant hepatocyte','Hepatocyte','Malignant hepatocyte','Immune/CAF',
                    'Malignant hepatocyte','Immune','Endothelial/CAF','Endothelial/CAF','Malignant hepatocyte','Hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Immune/CAF','Hepatocyte')
P11T.celltype <-  c('Hepatocyte','Endothelial/CAF','Hepatocyte','Malignant hepatocyte','Malignant hepatocyte','Immune/CAF','Malignant hepatocyte','Endothelial/CAF','Immune/CAF','Hepatocyte','Hepatocyte',
                    'Malignant hepatocyte','Hepatocyte','Malignant hepatocyte','Immune/CAF','Hepatocyte','Hepatocyte','Hepatocyte','Immune/CAF','Malignant hepatocyte','Hepatocyte','Immune/CAF')

### Annotate the cell types to ST data
Annotate_celltype <- function(data, celltype){
  data <- SetIdent(data, value = data@meta.data$SCT_snn_res.1.5)
  data@meta.data$CAFs <- data@meta.data$CAF
  data@meta.data$CAF <- NULL
  names(celltype) <- levels(data)
  data <- RenameIdents(data, celltype)
  return(data)
}

P1T <- Annotate_celltype(P1T, P1T.celltype)
P3T <- Annotate_celltype(P3T, P3T.celltype)
P5T <- Annotate_celltype(P5T, P5T.celltype)
P7T <- Annotate_celltype(P7T, P7T.celltype)
P8T <- Annotate_celltype(P8T, P8T.celltype)
P9T <- Annotate_celltype(P9T, P9T.celltype)
P10T <- Annotate_celltype(P10T, P10T.celltype)
P11T <- Annotate_celltype(P1T1, P11T.celltype)

### Later, we used BayesSpace to obtain higher resolution of the ST data
### Save file
save.image(paste0(work_path, "Human_ST_annotated.RData"))
