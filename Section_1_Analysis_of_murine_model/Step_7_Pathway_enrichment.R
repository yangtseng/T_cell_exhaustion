###########################################
### Section 1, Anlaysis of murine model ###
###########################################

###################################
### Pathway enrichment analysis ###
###################################

set.seed(1234)
work_path = "./"
source("requirements.R")

HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_TF6.rds"))

###################################################################
### Step 1, AUC score calculation of Hallmark and KEGG gene set ###
###################################################################

### We used AUCell to calculate the enrichment score of pathways for each cell
### We first built the rank via the gene expression profile
CellRank <- AUCell_buildRankings(as.matrix(HCC>tcell@assays$RNA@data))

### Get mouse gene sets data
### Here, we used KEGG and Hallmark gene sets
gs.go.h <- getGeneSets(species = "Mus musculus", library = 'H')
gs.go.kegg <- getGeneSets(species = 'Mus musculus', library = 'C2', subcategory = 'KEGG')

### Calculate the gene set AUC score for Hallmark gene set
cells_AUC_h <- AUCell_calcAUC(gs.go.h, CellRank, aucMaxRank=nrow(CellRank)*0.05)
save(cells_AUC_h, file = paste0(work_path, 'cells_AUC_h.RData'))

### Calculate the gene set AUC score for KEGG gene set
cells_AUC_kegg <- AUCell_calcAUC(gs.go.kegg, CellRank, aucMaxRank=nrow(CellRank)*0.05)
save(cells_AUC_kegg, file = paste0(work_path, 'cells_AUC_kegg.RData'))

################################################
### Step 2, Exploring the enrichment results ###
################################################

### We first added the enrichment results to seurat object
HCC.tcell[['KEGG']] <- CreateAssayObject(cells_AUC_kegg@assays@data@listData[["AUC"]])
HCC.tcell[['Hallmark']] <- CreateAssayObject(cells_AUC_h@assays@data@listData[["AUC"]])

### Next, we identify activated pathways in specific exhausted T cell clusters (cluster 1 and 4) compared to effector T cells
