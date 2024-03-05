###########################################
### Section 1, Anlaysis of murine model ###
###########################################

#######################################
### Transcription factor prediction ###
#######################################

set.seed(1234)
work_path = "./"
source("requirements.R")

load("murine_tcell_modulescore4.rds")

### We mainly followed the tutorial of SCENIC 
### http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

#############################################
### Step 1, Initialize settings of SCENIC ###
#############################################
org <- "mgi" # or hgnc, or dmel
dbDir <- paste0(work_path, "cisTarget-databases") # RcisTarget databases location
myDatasetTitle <- "TF_in_murine_model" # choose a name for your analysis

dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=12)
### save the initialize settings
saveRDS(scenicOptions, file=paste0(work_path, "int/scenicOptions.Rds"))

###################################################
### Step 2, Extraction of raw expression matrix ###
###################################################

### Raw counts from seurat object
exprMat <- as.matrix(HCC_tcell@assays[["RNA"]]@counts)
cellInfo <- data.frame(seuratClusters=Idents(HCC_tcell))
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
clusterInfo <- clusterInfo@cell.embeddings

