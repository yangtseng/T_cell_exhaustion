###################################################
### Section 2, Anlaysis of cell line like model ###
###################################################

######################################
### Module score and gene activity ###
######################################

set.seed(1234)

work_path = "./"
source(paste0(work_path, "requirements.R"))

### Here, we loaded the integrated multiome seurat object
tcell <- readRDS("Cellline_multiome1.rds")

##########################################
### Step 1, Calculate the module score ###
##########################################

### Set default assay as SCT
DefaultAssay(tcell) <- "SCT"

### Calculate the module score of exhausted and effectory gene set

