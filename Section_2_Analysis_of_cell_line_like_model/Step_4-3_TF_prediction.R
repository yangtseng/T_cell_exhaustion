###################################################
### Section 2, Anlaysis of cell line like model ###
###################################################

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


######################################
### Step 1, Import GRNboost output ###
######################################

GRNBoost_output <- read.delim(paste0(work_path, "int/1.1_grn_output.tsv"), header=FALSE)
colnames(GRNBoost_output) <- c("TF","Target","weight")
saveRDS(GRNBoost_output, file=paste0(work_path, "int/1.4_GRNboost_linkList.rds"))

#########################
### Step 2, RunSCENIC ###
#########################

### runSCENIC setting
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)

##################################################
### Step 3, Exploring SCENIC prediction result ###
##################################################

### Import regulon AUC data
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

### Extract the AUC score
regulonAUC_matrix <- regulonAUC@assays@data@listData[["AUC"]]
colnames(regulonAUC_matrix) <- regulonAUC@colData@rownames
rownames(regulonAUC_matrix) <- regulonAUC@NAMES

### load seurat object
tcell <- readRDS(paste0(work_path, "Cellline_multiome2.rds"))

### Create TF assay in seurat object
tcell[['TF']] <- CreateAssayObject(counts = regulonAUC_matrix)

############################################################
### Step 4, Differential expression analysis with seurat ###
############################################################

### Set default assay to TF
DefaultAssay(tcell) <- 'TF'

### Scale the regulon AUC score
tcell <- ScaleData(tcell)

### Differentially expressed TF
TF.markers <- FindAllMarkers(tcell, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
Tex.TF.markers <- FindMarkers(tcell, only.pos = TRUE, ident.1 = c('Exhausted_T_72hr','Exhausted_T_96hr'), ident.2 = c('Active_T'), min.pct = 0.1, logfc.threshold = 0)

### Heatmap visualization [Supp. Fig. 7]
DoHeatmap(tcell, features = TF.amrkers$gene, slot = 'scale.data', raster = F) #+ scale_fill_viridis()

### save seurat object
write.csv(Tex.TF.markers, paste0(work_path, "Cellline_TF_markers.csv"))
saveRDS(tcell, paste0(work_path, "Cellline_TF3.rds"))
