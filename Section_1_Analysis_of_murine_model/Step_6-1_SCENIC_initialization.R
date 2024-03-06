###########################################
### Section 1, Anlaysis of murine model ###
###########################################

#######################################
### Transcription factor prediction ###
#######################################

set.seed(1234)
work_path = "./"
source("requirements.R")

HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_modulescore4.rds"))

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
exprMat <- as.matrix(HCC.tcell@assays[["RNA"]]@counts)

### Cell information
cellInfo <- data.frame(seuratClusters=Idents(HCC.tcell))
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo

### Extract cluster information from seurat object
clusterInfo <- HCC.tcell@reductions[["umap"]]@cell.embeddings

###################################################
### Step 3, Construct the co-expression network ###
###################################################

### Genes overlap with RcisTarget
genesKept <- geneFiltering(exprMat, scenicOptions)

### There were 6905 genes available in RcisTarget
exprMat_filtered <- exprMat[genesKept, ]

### Construct the co-expression network
runCorrelation(exprMat_filtered, scenicOptions)

### log normalization
exprMat_filtered_log <- log2(exprMat_filtered + 1)

### Export the data for later prediction in R
saveRDS(exprMat_filtered_log, paste0(work_path, "int/exprMat_filtered_log.rds"))

### Export the data to use GRNboost in python
exportsForArboreto(exprMat_filtered_log, scenicOptions, dir = "int")
### We conducted the GRNboost on python, please refer to Section_1_Analysis_of_murine_model/Step_6-2_GRNboost.py
