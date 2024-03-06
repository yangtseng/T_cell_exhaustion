###################################################
### Section 2, Anlaysis of cell line like model ###
###################################################

########################################################################
### Module score, gene activity and differential expression analysis ###
########################################################################

set.seed(1234)

work_path = "./Cellline/"
source(paste0(work_path, "requirements.R"))

### Here, we loaded the integrated multiome seurat object
tcell <- readRDS("Cellline_multiome1.rds")

##########################################
### Step 1, Calculate the module score ###
##########################################

### Set default assay as SCT
DefaultAssay(tcell) <- "SCT"

### Calculate the module score of exhausted and effectory gene set

### T cell effectory gene set
eff <- list(c('Ifng','Ccl3','Ccl4','Prf1','Nkg7','Gzmb','Gzmk'))
tcell <- AddModuleScore(tcell, features = eff, name = 'effectory')

### T cell exhaustion gene set
exhaust <- list(c('Pdcd1','Cd244a', 'Lag3','Tigit','Eomes','Tox'))
tcell <- AddModuleScore(tcell, features = exhaust, name = 'exhaustion')

###########################################
### Step 2, Calculate the gene activity ###
###########################################

### Set default assay as peaks 
DefaultAssay(tcell) <- "peaks"

### Gene activity calculation
gene.activities <- GeneActivity(tcell)

### Create a new assay for gene activity storage
tcell[["GA"]] <- CreateAssayObject(counts = gene.activities)

### Normalization and scaling
DefaultAssay(exvivo240202) <- "GA"
tcell <- NormalizeData(tcell)
tcell <- ScaleData(tcell)

################################################
### Step 3, Differential expression analysis ###
################################################

### Identifying the differentially expressed genes (DEGs)
tcell_DEGs <- FindAllMarkers(tcell, assay = 'SCT', logfc.threshold = 0.4, only.pos = T)

### Identifying the differentially activated genes (based on gene activity)
tcell_DAGs <- FindAllMarkers(tcell, assay = 'GA', logfc.threshold = 0.3, only.pos = T)

### save files
write.csv(tcell_DEGs, paste0(work_path, "cell_line_like_DEGs.csv"))
write.csv(tcell_DAGs, paste0(work_path, "cell_line_like_DAGs.csv"))

saveRDS(tcell, paste0(work_path, "Cellline_multiome2.rds")
